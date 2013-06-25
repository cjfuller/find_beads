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
require 'ostruct'
require 'csv'

##
# Functions for segmenting and quantifying beads.
#
module FindBeads
  java_import Java::edu.stanford.cfuller.imageanalysistools.image.ImageFactory
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

    (next_dist_proj - closest_dist_proj < cutoff)
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
  # @param [OpenStruct, #centroids=] centroid_storage an optional object
  #   to which the bead centroids will be fed using #centroids=
  #   default=nil, causes them not to be set
  #
  def self.mask_from_image(im, opts, centroid_storage = nil)
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

    remapped_cens = {}
    cen_coord = ImageCoordinate[0,0,0,0,0]

    cens.each do |k, cen|
      cen_coord[:x] = cen[0]
      cen_coord[:y] = cen[1]
      val = final_mask[cen_coord]
      
      if val > 0 then
        remapped_cens[val] = cen
      end
    end

    cen_coord.recycle

    if centroid_storage then
      centroid_storage.centroids= remapped_cens
    end

    final_mask

  end

  ##
  # Swaps the x and y coordinates of an ImageCoordinate in place relative to a supplied origin.
  #
  # @param [ImageCoordinate] coord the coordinate whose components will be swapped.
  # @param [Numeric] cen_x the x-coordinate of the origin relative to which the swap will be calculated
  # @param [Numeric] cen_y the y-coordinate of the origin relative to which the swap will be calculated
  #
  # @return [ImageCoordinate] coord
  #
  def self.mirror_coord!(coord, cen_x, cen_y)
    rel_x = coord[:x] - cen_x
    rel_y = coord[:y] - cen_y
    coord[:x] = (rel_y + cen_x).round.to_i
    coord[:y] = (rel_x + cen_y).round.to_i
    coord
  end

  ##
  # Computes the two-channel correlation for a single bead.  
  #
  # @param [Image] mask the mask of labeled, segmented beads
  # @param [Image] im the image of the beads (single z-section, multiple channels)
  # @param [Integer] ch0 the (0-based) index of the first channel of the correlation
  # @param [Integer] ch1 the (0-based) index of the second channel of the correlation
  # @param [Array] cen an array containing the x,y coordinates of the centroid of the current bead
  # @param [Integer] id the label of the current bead in the mask
  # @param [Numeric] bead_radius the radius of a bead
  # @param [Boolean] do_normalization whether to normalize to the geometric mean of the autocorrelations of the two channels
  #
  # @return [Hash] a hash containing components :norm_corr (the background-subtracted, normalized correlation), :corr (the
  #                non-normalized, non-subtracted two-channel correlation), and :bg_corr (the non-normalized background 
  #                two-channel correlation, calculated by rotating one channel of the bead 90 degrees around the bead centroid)
  #
  def self.compute_single_bead_correlation(mask, im, ch0, ch1, cen, id, bead_radius, do_normalization)
    normalization = if do_normalization then
                      auto_00 = compute_single_bead_correlation(mask, im, ch0, ch0, cen, id, bead_radius, false)[:norm_corr]
                      auto_11 = compute_single_bead_correlation(mask, im, ch1, ch1, cen, id, bead_radius, false)[:norm_corr]
                      begin
                        Math.sqrt(auto_00 * auto_11)*(auto_00.abs/auto_00)
                      rescue Math::DomainError
                        Float::NAN
                      end
                    else
                      1
                    end

    box_lower = ImageCoordinate[cen[0] - bead_radius - 1, cen[1] - bead_radius - 1,0,0,0]
    box_upper = ImageCoordinate[cen[0] + bead_radius + 2, cen[1] + bead_radius + 2,1,1,1]
    mask.setBoxOfInterest(box_lower, box_upper)
    
    ch_coord = ImageCoordinate[0,0,0,0,0]
    corr_sum = 0
    bg_sum = 0
    count = 0
    bg_count = 0

    mask.each do |ic|
      next unless mask[ic] == id
      ch_coord.setCoord(ic)
      ch_coord[:c] = ch0
      val = im[ch_coord]
      ch_coord[:c] = ch1
      corr_sum += val * im[ch_coord]
      count += 1
      mirror_coord!(ch_coord, cen[0], cen[1])
      if im.inBounds(ch_coord) then
        bg_sum += val * im[ch_coord]
        bg_count += 1
      end
    end
      
    mask.clearBoxOfInterest
    box_lower.recycle
    box_upper.recycle
    ch_coord.recycle

    {norm_corr: (corr_sum/count - bg_sum/bg_count)/normalization , corr: corr_sum/count, bg_corr: bg_sum/bg_count}

  end

  ##
  # Runs a correlation analysis between two channels of the image.  This is done by directly
  # computing the correlation for each bead and then computing an effective background by
  # swapping the x- and y- axes of one channel within each bead.
  #
  # @param [Image] mask the mask of the beads
  # @param [Integer] ch0 the first channel for the correlation
  # @param [Integer] ch1 the second channel for the correlation
  # @param [Hash] cens a hash mapping bead labels in the mask to the centroids of the beads
  # @param [Hash] opts the options hash
  #
  # @return [Hash] a hash mapping :norm_corr, :corr, and :bg_corr to a hash mapping bead labels to their measurements
  #
  def self.compute_correlation(mask, im, ch0, ch1, cens, opts)
    result = {norm_corr: {}, corr: {}, bg_corr: {}}
    cens.each do |id, cen|
      corr_for_bead = compute_single_bead_correlation(mask, im, ch0, ch1, cen, id, opts[:beadradius], true)
      result[:norm_corr][id] = corr_for_bead[:norm_corr]
      result[:corr][id] = corr_for_bead[:corr]
      result[:bg_corr][id] = corr_for_bead[:bg_corr]
    end
    result
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
  # Organizes the correlation data into a csv string where each row is a bead and each column is a measurement.
  # Also creates a header row.
  #
  # @param [Hash] corr_output a hash formatted like the output of FindBeads.compute_correlation
  # 
  # @return [String] a string containing the same data in csv format.
  #
  def self.format_correlation_output(corr_output)
    header_row = ["bead index", "normalized, corrected correlation", "correlation", "background correlation"]
    
    CSV.generate do |csv|
      csv << header_row
      
      corr_output[:corr].each_key do |id|
        csv << [id, corr_output[:norm_corr][id], corr_output[:corr][id], corr_output[:bg_corr][id]]
      end
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

    cens = if opts and opts[:correlation_channels] then
             OpenStruct.new
           else
             nil
           end

    im = RImageAnalysisTools.get_image(fn)
    mask = mask_from_image(im, opts, cens)
    proj = Java::edu.stanford.cfuller.imageanalysistools.frontend.MaximumIntensityProjection.projectImage(im)
    ims = proj.splitChannels
    is = ImageSet.new(ParameterDictionary.emptyDictionary)

    ims.each do |imc|
      is.addImageWithImage(imc)
    end

    outdat = if opts and opts[:correlation_channels] then
               q = compute_correlation(mask, proj, opts[:correlation_channels][0], opts[:correlation_channels][1], cens.centroids, opts)
               format_correlation_output(q)
             else
               met = IntensityPerPixelMetric.new
               q = met.quantify(mask, is)
               Java::edu.stanford.cfuller.imageanalysistools.frontend.LocalAnalysis.generateDataOutputString(q, nil).strip.gsub(" ", ",")
             end

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
      opt :max_threads, "Maximum number of paralell execution threads", :type => :integer, :default => DEFAULT_THREADS
      opt :beadradius, "Radius of the bead in pixels", :type => :float, :default => DEFAULT_BEAD_RADIUS
      opt :correlation_channels, "Runs correlation between the specified two comma separated channels", :type => :string, :default => nil
    end

    if opts[:correlation_channels] then
      opts[:correlation_channels] = opts[:correlation_channels].gsub("\s+", "").split(",").map(&:to_i)
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

