require 'lib/find_beads/version'

Gem::Specification.new do |g|

  g.name = "find_beads"
  g.version = FindBeads::VERSION
  g.date = '2013-02-02'
  g.summary = "Segments and quantifies beads in microscopy images"  
  g.description = g.summary
  g.authors = ['Colin J. Fuller']
  g.email = 'cjfuller@gmail.com'
  g.homepage = 'http://github.com/cjfuller/find_beads'
  g.add_runtime_dependency 'rimageanalysistools'
  g.add_runtime_dependency 'trollop'
  g.files = Dir['lib/**/*.rb', 'bin/**/*']
  g.executables << 'find_beads'
  g.license = 'MIT'
  g.platform = 'java'
  g.requirements = 'jruby'

end
  
