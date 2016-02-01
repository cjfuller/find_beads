import java.io.File
import scala.collection.JavaConverters._
import scala.concurrent._
import ExecutionContext.Implicits.global

import scopt.OptionParser

import io.cjf.findbeads.FindBeads
import io.cjf.findbeads.config._

object Main {
  def main(args: Array[String]) {
    val config = handleArgs(args)

    config match {
      case Some(c) => {
        c.toProcess match {
          case None => require(false, "Must provide either a file or directory")
          case Some(SingleFile(f)) => FindBeads.processFile(f, c)
          case Some(Dir(d)) => {
            require(d.exists(), "Directory does not exist.")
            d.listFiles().filter(_.isFile).filter(!_.getName().startsWith(".")).map(
              (f) => Future {FindBeads.processFile(f, c)}
            ).foreach((f) => Await.result(f, duration.Duration.Inf))
          }
        }
      }
      case None => None
    }
  }

  def handleArgs(args: Array[String]): Option[Config] = {
    val parser = new scopt.OptionParser[Config]("find_beads") {
      head("find beads v1")
      opt[File]('f', "file") action { (x, c) => c.copy(toProcess = Some(SingleFile(x)))
      } text("file to process")
      opt[File]('d', "dir") action { (x, c) => c.copy(toProcess = Some(Dir(x)))
      } text("directory to process")
      opt[Int]("segchannel") action { (x, c) => c.copy(segchannel = x)
      } text(s"Channel on which to segment (0-indexed); default: ${Const.DEFAULT_SEG_CH}")
      opt[Int]("segplane") action { (x, c) => c.copy(segplane = x)
      } text(s"Plane on which to segment (0-indexed); default: ${Const.DEFAULT_SEG_PL}")
      opt[Double]("beadradius") action { (x, c) => c.copy(beadradius = x)
      } text(s"Radius of a bead in pixels; default: ${Const.DEFAULT_BEAD_RADIUS}")
      help("help") text("prints this usage text")
    }

    parser.parse(args, Config())
  }
}
