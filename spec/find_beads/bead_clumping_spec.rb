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

require 'find_beads/bead_clumps'


describe FindBeads::BeadClumping do

	before :each do 

		@points = [[0, 0], [1, 0], [2, 0], [0, -1], [1, -1], #group 1
					[0, 3], [1, 3], [1, 4], #group 2
					[-2, 0], [-2, 1], [-2, 2], [-3, 1], [-3, 2]] #group 3

		@cutoff = 1.3

	end

	it "should be able to find the labels of beads within a certain cutoff distance" do

		component = FindBeads::BeadClumping.get_points_within_distance([0, 3], @points, @cutoff)

		component.include?([1, 3]).should be_true
		component.include?([1, 4]).should be_false

	end

	it "should be able to construct a graph of interconnected beads within a certain cutoff distance" do

		g = FindBeads::BeadClumping.construct_adjacency_graph(@points, @cutoff)

		g.should be_connected(1, 2)

		g.should_not be_connected(1, 6)

	end

	it "should be able to calculate the number of disjoint subgraphs" do

		g = FindBeads::BeadClumping.construct_adjacency_graph(@points, @cutoff)

		FindBeads::BeadClumping.count_disjoint_subgraphs(g).should eq(3)

	end

end
