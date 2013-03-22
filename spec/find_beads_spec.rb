#--
# Copyright (c) 2013 Colin J. Fuller
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the Software), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#++

require 'spec_helper'

require 'rimageanalysistools'
require 'rimageanalysistools/get_image'
require 'rimageanalysistools/image_shortcuts'

require 'find_beads'


describe FindBeads do

	describe "image processing operations" do

		before :all do

			@im = RImageAnalysisTools.get_image(IMG_FN)
			@im_out = RImageAnalysisTools.get_image(IMG_OUT_FN)

		end

		it "should correctly segment a bead image", :slow do

			FindBeads.mask_from_image(@im, test_cmd_line_args).should look_like(@im_out)

		end

		it "should correctly calculate the centroids of objects in a mask" do

			im_test = Java::edu.stanford.cfuller.imageanalysistools.image.ImageFactory.create(ImageCoordinate[11,11,1,1,1], 1.0)

			expected = { 1.0 => [5.0, 5.0] }

			FindBeads.centroids(im_test).should eq expected


			expected_realimage = {1.0=>[97.11400110680687, 44.54676258992806], 2.0=>[49.38269550748752, 47.09317803660566], 3.0=>[468.9717673630717, 247.84246188594014], 4.0=>[565.7143668370244, 261.541737649063], 5.0=>[508.3385689354276, 268.6038394415358], 6.0=>[628.4727168294083, 282.6105686387134], 7.0=>[710.3073367995379, 298.41247833622185], 8.0=>[538.1956652618904, 299.48886213124626], 9.0=>[582.6725768321513, 303.99881796690306], 10.0=>[666.8053613053613, 306.8706293706294], 11.0=>[747.9779536461277, 322.5166760881854], 12.0=>[619.7359413202934, 326.7658924205379], 13.0=>[543.5572082379863, 343.6567505720824], 14.0=>[655.6335439402642, 351.8012636415853], 15.0=>[587.2113402061856, 357.0171821305842], 16.0=>[901.771270718232, 466.1265193370166], 17.0=>[298.52459922609177, 495.3819789939193], 18.0=>[287.70550724637684, 575.3298550724637], 19.0=>[250.37528735632185, 599.5137931034483], 20.0=>[324.5933294324166, 599.8665886483324], 21.0=>[398.14368147720717, 623.6093479515291], 22.0=>[354.07134292565945, 633.3489208633093], 23.0=>[226.67504283266706, 638.0868075385494], 24.0=>[310.8706438275251, 643.9031305375074], 25.0=>[269.30094228504123, 654.228504122497], 26.0=>[384.36026331538, 666.3734290843806], 27.0=>[343.23244837758114, 683.0259587020649], 28.0=>[238.8558352402746, 688.0], 29.0=>[305.55523255813955, 706.0691860465116], 30.0=>[406.57012542759406, 706.2377423033067], 31.0=>[451.10597519729424, 714.4909808342728], 32.0=>[263.8544194107452, 725.6891969959561], 33.0=>[348.44707207207205, 728.5557432432432], 34.0=>[300.5765969474279, 752.9745618993782]}

			FindBeads.centroids(@im_out).should eq expected_realimage

		end

	end

	it "should correctly perform projection of a point onto a vector" do

		origin = [0.0, 0.0]
		point_on_line = [3.0, 0.0]
		point_to_project = [7.0, 2.0]

		expected = [7.0, 0.0]

		FindBeads.project_point_onto_vector(origin, point_to_project, point_on_line).should eq expected

	end

	it "should correctly calculate region borders of a voronoi diagram" do 

		points = [[-2, 2], [2, 2]]
		ic_on = ImageCoordinate[0,2,0,0,0]
		ic_off = ImageCoordinate[1,2,0,0,0]

		FindBeads.is_on_voronoi_border?(points, ic_off).should be_false
		FindBeads.is_on_voronoi_border?(points, ic_on).should be_true

	end

	it "should calculate a reasonable maximum bead size for a given radius" do

		rad = 10

		max_size = FindBeads.calculate_max_size_from_radius(rad)

		max_size.should be < ((2*rad + 1)**2)
		max_size.should be > (Math::PI * rad**2)

	end

	it "should calculate a reasonable minimum bead size for a given radius" do

		rad = 10

		min_size = FindBeads.calculate_min_size_from_radius(rad)

		min_size.should be < (Math::PI * rad**2)
		min_size.should be > (0.25 * Math::PI * rad**2)

	end

end
