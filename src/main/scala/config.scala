package io.cjf.findbeads.config

import java.io.File

object Const {
  val DEFAULT_SEG_CH = 2
  val DEFAULT_SEG_PL = 8
  val DEFAULT_BEAD_RADIUS = 24.0
}

sealed abstract class FileOrDir
case class SingleFile(f: File) extends FileOrDir
case class Dir(f: File) extends FileOrDir

case class Config(
  toProcess: Option[FileOrDir] = None,
  segchannel: Int = Const.DEFAULT_SEG_CH,
  segplane: Int = Const.DEFAULT_SEG_PL,
  beadradius: Double = Const.DEFAULT_BEAD_RADIUS
)
