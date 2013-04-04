#--
# /* ***** BEGIN LICENSE BLOCK *****
#  * 
#  * Copyright (c) 2013 Colin J. Fuller
#  * 
#  * Permission is hereby granted, free of charge, to any person obtaining a copy
#  * of this software and associated documentation files (the Software), to deal
#  * in the Software without restriction, including without limitation the rights
#  * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  * copies of the Software, and to permit persons to whom the Software is
#  * furnished to do so, subject to the following conditions:
#  * 
#  * The above copyright notice and this permission notice shall be included in
#  * all copies or substantial portions of the Software.
#  * 
#  * THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  * SOFTWARE.
#  * 
#  * ***** END LICENSE BLOCK ***** */
#++

require 'rimageanalysistools'
require 'rimageanalysistools/get_image'
require 'rimageanalysistools/image_shortcuts'
require 'rimageanalysistools/create_parameters'
require 'rimageanalysistools/graythresh'

require 'trollop'

require 'matrix'

java_import Java::edu.stanford.cfuller.imageanalysistools.image.ImageFactory

##
# Functions for segmenting and quantifying beads.
#
module FindBeads


  java_import Java::edu.stanford.cfuller.imageanalysistools.image.ImageCoordinate
  java_import Java::edu.stanford.cfuller.imageanalysistools.image.ImageSet
  java_import Java::edu.stanford.cfuller.imageanalysistools.meta.parameters.ParameterDictionary
  java_import Java::edu.stanford.cfuller.imageanalysistools.filter.MaximumSeparabilityThresholdingFilter
  java_import Java::edu.stanford.cfuller.imageanalysistools.filter.LabelFilter
  java_import Java::edu.stanford.cfuller.imageanalysistools.filter.SizeAbsoluteFilter
  java_import Java::edu.stanford.cfuller.imageanalysistools.filter.VoronoiFilter
  java_import Java::edu.stanford.cfuller.imageanalysistools.filter.MaskFilter
  java_import Java::edu.stanford.cfuller.imageanalysistools.metric.IntensityPerPixelMetric
  java_import Java::edu.stanford.cfuller.imageanalysistools.image.Histogram

  DEFAULT_SEG_CH = 2
  DEFAULT_SEG_PL = 8
  DEFAULT_BEAD_RADIUS = 24.0
  DEFAULT_THREADS = 1

  ##
  # Finds the centroid of each unique-greylevel region in a mask.
  #
  # @param [Image] mask  the mask in which the regions appear.  0 denotes background and will not be counted.
  # @return [Hash] a hash where keys are the greylevels of each region in the mask, and the values are
  #  two-element arrays containing the x,y-coordinates of the centroids of these regions.
  #
  def self.centroids(mask)

    cens = {}

    mask.each do |ic|

      next unless mask[ic] > 0
      
      cens[mask[ic]] = [0.0, 0.0] unless cens[mask[ic]]

      cens[mask[ic]][0] += ic[:x]
      cens[mask[ic]][1] += ic[:y]
      
    end

    h = Histogram.new(mask)

    cens.each_key do |k|

      cens[k].map! { |e| e / h.getCounts(k) }

    end

    cens

  end


  ##
  # Projects a point in space onto a specified line.
  #
  # @param [Array] origin the point in space the will serve as the origin 
  #  for the purposes of the projection (this should be on the line)
  # @param [Array] point_to_project the point in space that will be projected
  #  this should be in the same coordinate system in which the origin is
  #  specified, not relative to the origin
  # @param [Array] point_on_line another point on the line specified in the
  #  same coordinate system in which the origin is specified, not relative to
  #  the origin
  #
  # @return [Array] the projected point (in the same coordinate system as
  #  the other points were specified)
  #
  def self.project_point_onto_vector(origin, point_to_project, point_on_line)

    unit_vec = (Vector[*point_on_line] - Vector[*origin]).normalize

    proj_vec = Vector[*point_to_project] - Vector[*origin]

    projected = unit_vec * (proj_vec.inner_product(unit_vec)) + Vector[*origin]

    projected.to_a

  end


  ##
  # Checks if a given coordinate would be approximately on the boundary between two regions of a 
  # Voronoi diagram of constructed from a set of points.  The approximation is calculated such that no
  # two regions in the Voronoi diagram would be 8-connected.
  #
  # @param [Array] points an array of two-element arrays containing the x,y coordinates of the points
  #  on which the diagram is calculated.
  # @param [ImageCoordinate] ic an ImageCoordinate specifying the location to check.
  # @return [Boolean] whether the point is on the border between two regions.
  #
  def self.is_on_voronoi_border?(points, ic)

    x = ic[:x]
    y = ic[:y]

    closest_index = 0
    next_index = 0
    closest_dist = Float::MAX
    next_dist = Float::MAX
    

    points.each_with_index do |p, i|

      dist = Math.hypot(p[0] - x, p[1] - y)

      if dist < closest_dist then

        next_dist = closest_dist

        next_index = closest_index

        closest_dist = dist

        closest_index = i

      elsif dist < next_dist then

        next_dist = dist

        next_index = i

      end

    end

    proj_point = project_point_onto_vector(points[closest_index], [x,y], points[next_index])

    next_dist_proj = Math.hypot(points[next_index][0]-proj_point[0], points[next_index][1]-proj_point[1])
    closest_dist_proj = Math.hypot(points[closest_index][0]-proj_point[0], points[closest_index][1]-proj_point[1])


    cutoff = 1.01*Math.sqrt(2)

    if next_dist_proj - closest_dist_proj < cutoff then

      true

    else 

      false

    end

  end

  ##
  # Recursively thresholds regions in a supplied mask and image using the method described in 
  # Xiong et al. (DOI: 10.1109/ICIP.2006.312365).
  #
  # @param [ParameterDictionary] p a ParameterDictionary specifying max_size and min_size parameters,
  #  which control the maximum size of regions before they are recursively thresholded to break them up,
  #  and the minimum size of regions before they are discarded.
  # @param [Image] im the original image being segmented.  This will not be modified.
  # @param [Image] mask the initial segmentation mask for the supplied image.  Regions in this mask may
  #  be divided up or discarded.
  # @return [void]
  # 
  def self.recursive_threshold(p, im, mask)

    h = Histogram.new(mask)

    changed = false

    discard_list = {}

    1.upto(h.getMaxValue) do |i|

      if h.getCounts(i) > p[:max_size].to_i then

        values = []

        im.each do |ic|

          if mask[ic] == i then

            values << im[ic]

          end

        end

        thresh = RImageAnalysisTools.graythresh(values)

        im.each do |ic|

          if mask[ic] == i and im[ic] <= thresh then

            mask[ic] = 0

            changed = true

          end

        end

      elsif h.getCounts(i) > 0 and h.getCounts(i) < p[:min_size].to_i then

        discard_list[i] = true

      end

    end

    im.each do |ic|

      if discard_list[im[ic]] then

        mask[ic] = 0

        changed = true
        
      end

    end


    if changed then

      lf = LabelFilter.new

      lf.apply(mask)

      recursive_threshold(p, im, mask)

    end

  end

  ##
  # Caclulates the maximum allowed size of a bead from the supplied radius.
  # This is set to be slightly larger than a circle of that radius.
  #
  # @param [Fixnum] rad the radius of the bead in units of pixels.
  #
  def self.calculate_max_size_from_radius(rad)

    ((rad+1)**2 * 3.2).to_i

  end

  ##
  # Calculates the minimum allowed size of a bead from the supplied radius.
  # This is set to be slightly smaller than a third of a circle of that radius.
  # (Note that this is smaller than any of the returned regions should be, but making
  # the cutoff this small is useful for dividing up clumps of beads where several rounds
  # or recursive thresholding may make the regions quite small temporarily.)
  #
  # @param @see #calculate_max_size_from_radius
  #
  def self.calculate_min_size_from_radius(rad)

    (0.96* (rad+1)**2).to_i

  end

  ##
  # Generates a segmented mask of beads from an image.
  #
  # @param [Image] im the image to segment
  # @param [Hash] opts a hash of commandline arguments.
  #
  def self.mask_from_image(im, opts)

    seg_ch = nil
    seg_pl = nil
    rad = nil

    if opts then
      seg_ch = opts[:segchannel]
      seg_pl = opts[:segplane]
      rad = opts[:beadradius]
    else
      seg_ch = DEFAULT_SEG_CH
      seg_pl = DEFAULT_SEG_PL
      rad = DEFAULT_BEAD_RADIUS
    end
    
    min_size = calculate_min_size_from_radius(rad)
    max_size = calculate_max_size_from_radius(rad)
    
    sizes = ImageCoordinate.cloneCoord(im.getDimensionSizes)

    sizes[:c] = 1
    sizes[:z] = 1

    im0 = ImageCoordinate.createCoordXYZCT(0,0,0,0,0)

    im0[:c] = seg_ch
    im0[:z] = seg_pl

    to_seg = im.subImage(sizes, im0).writableInstance

    p = RImageAnalysisTools.create_parameter_dictionary(min_size: min_size, max_size: max_size)

    im_cp = ImageFactory.create_writable(to_seg)

    mstf = MaximumSeparabilityThresholdingFilter.new

    lf = LabelFilter.new

    saf = SizeAbsoluteFilter.new

    filters = [] 

    filters << mstf

    filters << lf

    filters.each do |f|

      f.setParameters(p)
      f.setReferenceImage(im_cp)
      f.apply(to_seg)

    end

    recursive_threshold(p, im_cp, to_seg)

    saf.setParameters(p)
    
    saf.apply(to_seg)

    cens = centroids(to_seg)

    final_mask = ImageFactory.create_writable(to_seg)

    radius = rad

    final_mask.each do |ic|

      final_mask[ic] = 0
      
    end

    final_mask.each do |ic|

      x = ic[:x]
      y = ic[:y]

      cens.each_key do |k|

        if Math.hypot(cens[k][0] - x, cens[k][1] - y) <= radius then

          final_mask[ic] = k

        end

      end

    end


    final_mask.each do |ic|

      next unless final_mask[ic] > 0

      if is_on_voronoi_border?(cens.values, ic) then

        final_mask[ic] = 0

      end

    end

    lf.apply(final_mask)
    
    saf.apply(final_mask)

    lf.apply(final_mask)

    final_mask

  end

  ##
  # Writes the output data and mask to files.
  #
  # @param [String] fn_orig  the original filename of the image being segmented/quantified.
  # @param [String] quant_str the quantification data to write
  # @param [Image] mask the mask to write
  #
  def self.write_output(fn_orig, quant_str, mask)

    mask_dir = "output_mask"

    quant_dir = "quantification"

    mask_ext = "_mask.ome.tif"

    quant_ext = "_quant.txt"

    dir = File.dirname(fn_orig)

    base = File.basename(fn_orig)

    base = base.gsub(".ome.tif", "")

    mask_dir = File.expand_path(mask_dir, dir)

    quant_dir = File.expand_path(quant_dir, dir)

    Dir.mkdir(mask_dir) unless Dir.exist?(mask_dir)
    Dir.mkdir(quant_dir) unless Dir.exist?(quant_dir)

    mask.writeToFile(File.expand_path(base + mask_ext, mask_dir))

    File.open(File.expand_path(base + quant_ext, quant_dir), 'w') do |f|

      f.puts(quant_str)

    end

  end

  ##
  # Processes a single file, which consists of creating a mask, quantifying regions, and writing output.
  #
  # @param [String] fn the filename of the image to process
  # @param [Hash] opts a hash of command line options.
  #
  def self.process_file(fn, opts=nil)

    puts "processing #{fn}"

    im = RImageAnalysisTools.get_image(fn)
    
    mask = mask_from_image(im, opts)

    proj = Java::edu.stanford.cfuller.imageanalysistools.frontend.MaximumIntensityProjection.projectImage(im)

    ims = proj.splitChannels

    is = ImageSet.new(ParameterDictionary.emptyDictionary)

    ims.each do |imc|

      is.addImageWithImage(imc)

    end

    met = IntensityPerPixelMetric.new

    q = met.quantify(mask, is)

    outdat = Java::edu.stanford.cfuller.imageanalysistools.frontend.LocalAnalysis.generateDataOutputString(q, nil)

    write_output(fn, outdat, mask)

  end

  ##
  # Runs the bead finding on a file or directory, and grabs options from the command line.
  #
  def self.run_find_beads

    opts = Trollop::options do

      opt :dir, "Directory to process", :type => :string
      opt :file, "File to process", :type => :string
      opt :segchannel, "Channel on which to segment (0-indexed)", :type => :integer, :default => DEFAULT_SEG_CH
      opt :segplane, "Plane on which to segment (0-indexed)", :type => :integer, :default => DEFAULT_SEG_PL
      opt :max_threads, "Maximum number of parallel execution threads", :type => :integer, :default => DEFAULT_THREADS
      opt :beadradius, "Radius of the bead in pixels", :type => :float, :default => DEFAULT_BEAD_RADIUS

    end

    if opts[:dir] then

      fod = opts[:dir]

      sleep_time_s = 0.5
      threads = []

      Dir.foreach(fod) do |f|

        until threads.count { |t| t.alive? } < opts[:max_threads] do
          sleep sleep_time_s
        end

        fn = File.expand_path(f, fod)

        if File.file?(fn) then

          begin

            threads << Thread.new do 

              process_file(fn, opts)

            end

          rescue Exception => e

            puts "Unable to process #{fn}:"
            puts e.message

          end
          
        end

      end

      threads.each { |t| t.join }

    end



    if opts[:file] then

      process_file(opts[:file], opts)
      
    end

  end

end

