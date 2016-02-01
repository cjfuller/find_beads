#Introduction

Script to segment and quantify beads in microscopy images.

#Installation

You'll need java; this has only been tested on Java 8.

Download the executable jar file here: https://storage.googleapis.com/find-beads-dist/find-beads-1.0

On OSX / linux, you can run `chmod u+x find-beads-1.0` to run it as if it's an executable itself.

#Running
`find-beads-1.0 --help`
or
`java -jar find-beads-1.0 --help`

for usage instructions.

#Building from source
You'll need Scala (tested only with 2.11) and sbt.

`sbt assembly` will build an executable jar with all dependencies.

#License

Find_beads is distrubuted under the MIT/X11 license.  See the full text in the LICENSE file for more details.  Some dependencies are distributed under different licenses.  See their project pages for more details.



