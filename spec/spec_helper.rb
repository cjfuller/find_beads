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

require 'rimageanalysistools'
require 'rimageanalysistools/image/io/url_image_reader'

Java::org.apache.log4j.Logger.getRootLogger().setLevel(Java::org.apache.log4j.Level::OFF)

RSpec.configure do |c|
	c.treat_symbols_as_metadata_keys_with_true_values = true
end

Dir[File.dirname(__FILE__) + "/support/**/*.rb"].each { |f| require f }

IMG_FN = URI("https://s3-us-west-1.amazonaws.com/find-beads/test_bead_input.ome.tif")

IMG_OUT_FN = URI("https://s3-us-west-1.amazonaws.com/find-beads/test_bead_output.ome.tif")

def test_cmd_line_args

	{segchannel: 0, segplane: 4, beadradius: 24.0}

end
