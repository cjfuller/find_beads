import sbtassembly.AssemblyPlugin.defaultShellScript

name := "find-beads"

version := "1.0"
scalaVersion := "2.11.7"

libraryDependencies += "com.github.scopt" %% "scopt" % "3.3.0"
libraryDependencies += "org.scalatest" %% "scalatest" % "2.2.6" % "test"

assemblyOption in assembly := (assemblyOption in assembly).value.copy(prependShellScript = Some(defaultShellScript))
assemblyJarName in assembly := s"${name.value}-${version.value}"
