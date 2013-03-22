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

require 'set'

module FindBeads

	module BeadClumping

		class Edge

			def initialize(v0, v1)
				@v0 = v0
				@v1 = v1
			end

			attr_accessor :v0, :v1

			def connected_to?(v)
				v == @v0 or v == @v1
			end

			def ==(e)
				(e.v0 == @v0 and e.v1 == @v1) or (e.v1 == @v0 and e.v0 == @v1)
			end

			alias eql? ==

		end

		class Vertex

			def initialize(label)
				@edges = Set.new
				@label = label
			end

			def connect(v)
				e = Edge.new(self, v)
				@edges << e
			end

			def add_edge(e)
				@edges << e if (e.v0 == self or e.v1 == self)
			end

			def connected_to?(v)
				e.Edge.new(self, v)
				@edges.include?(e)
			end

			def include?(e)
				@edges.include?(e)
			end

			def each(&b)
				@edges.each &b
			end

			attr_accessor :edges, :label

		end

		class Graph

			def initialize
				@vertices = {}
			end

			def add_vertex(v)
				@vertices[v.label]= v
			end

			def new_vertex(l)
				@vertices[l] = Vertex.new(l)
			end

			def vertex(label)
				@vertices[label]
			end

			def each(&b)
				@vertices.each &b
			end

			alias [] vertex

			def direct_connected?(l0, l1)
				self[l0].connected_to? self[l1]
			end

			def add_edge_between(l0, l1)
				self[l0].connect(self[l1])
			end

			def connected?(l0, l1)
				sg = get_connected_subgraph(l0)
				sg.include?(self[l1])
			end

			def build_connected_dfs(current, visited)

				visited << current

				current.each do |e|

					[:v0, :v1].each do |vi|

						v = e.send(vi)

						if v.nil?

							puts e.inspect

							raise

						end

						unless visited.include?(v) then

							build_connected_dfs(v, visited)

						end

					end

				end

				visited

			end

			def get_connected_subgraph(l0)

				build_connected_dfs(self[l0], Set.new)

			end

		end


		def self.get_points_within_distance(point, points_list, distance_cutoff)

			close = []

			points_list.each do |p|

				next if point == p

				if (Vector[*point] - Vector[*p]).norm < distance_cutoff then

					close << p

				end

			end

			close

		end

		def self.construct_adjacency_graph(points, distance_cutoff)

			g = Graph.new

			points.each_with_index do |p, i|

				g.new_vertex(i+1)

			end

			points.each_with_index do |p, i|

				close = get_points_within_distance(p, points, distance_cutoff)

				close.each do |c|

					g.add_edge_between(i+1, points.index(c) + 1)

				end

			end

			g

		end


		def self.count_disjoint_subgraphs(graph)

			subgraphs = Set.new

			graph.each do |l, vert|

				subgraphs << graph.get_connected_subgraph(l)

			end

			subgraphs.size

		end

	end

end


