
require 'rimageanalysistools'
require 'rimageanalysistools/get_image'
require 'rimageanalysistools/image_shortcuts'
require 'rimageanalysistools/create_parameters'
require 'rimageanalysistools/graythresh'

require 'edu/stanford/cfuller/imageanalysistools/resources/common_methods'

require 'trollop'

include IATScripting

java_import Java::edu.stanford.cfuller.imageanalysistools.image.ImageCoordinate
java_import Java::edu.stanford.cfuller.imageanalysistools.image.ImageSet

java_import Java::edu.stanford.cfuller.imageanalysistools.meta.parameters.ParameterDictionary

java_import Java::edu.stanford.cfuller.imageanalysistools.filter.RecursiveMaximumSeparabilityFilter
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


def centroids(mask)

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

def is_on_voronoi_border?(points, ic)

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

    end

  end

  cutoff = 2*Math.sqrt(2)

  if next_dist - closest_dist < cutoff then

    true

  else 

    false

  end

end


def recursive_threshold(p, im, mask)

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

def calculate_max_size_from_radius(rad)

  ((rad+1)**2 * 3.2).to_i

end

def calculate_min_size_from_radius(rad)

  (0.96* (rad+1)**2).to_i

end

def mask_from_image(im, opts)

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

  mstf = MaximumSeparabilityThresholdingFilter.new

  im_cp = writable_image_copy(to_seg)

  rmsf = RecursiveMaximumSeparabilityFilter.new

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

  final_mask = writable_image_copy(to_seg)

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

def write_output(fn_orig, quant_str, mask)

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

  

def process_file(fn, opts=nil)

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

  

if __FILE__ == $0 then

  opts = Trollop::options do

    opt :dir, "Directory to process", :type => :string
    opt :file, "File to process", :type => :string
    opt :segchannel, "Channel on which to segment (0-indexed)", :type => :integer, :default => DEFAULT_SEG_CH
    opt :segplane, "Plane on which to segment (0-indexed)", :type => :integer, :default => DEFAULT_SEG_PL
    opt :beadradius, "Radius of the bead in pixels", :type => :float, :default => DEFAULT_BEAD_RADIUS

  end

  if opts[:dir] then

    fod = opts[:dir]

    Dir.foreach(fod) do |f|

      fn = File.expand_path(f, fod)

      if File.file?(fn) then

        begin

          process_file(fn, opts)

        rescue Exception => e

          puts "Unable to process #{fn}:"
          puts e.message

        end
          
      end

    end

  end

  if opts[:file] then

    process_file(opts[:file], opts)
    
  end

end
