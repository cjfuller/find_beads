package io.cjf.findbeads.test

import java.io.{File, FileOutputStream, BufferedOutputStream, BufferedInputStream}
import java.net.{URL, URLConnection}
import scala.collection.mutable.ArrayBuffer
import scala.collection.JavaConverters._

import edu.stanford.cfuller.imageanalysistools.image.{Image, ImageCoordinate, ImageFactory}
import edu.stanford.cfuller.imageanalysistools.image.io.ImageReader
import org.scalatest._
import matchers._

import io.cjf.findbeads.FindBeads
import io.cjf.findbeads.config.Config

object SlowTest extends Tag("io.cjf.findbeads.SlowTest")

trait LookLikeMatcher {
  class LookLikeMatcherCls(expectedImage: Image) extends Matcher[Image] {
    def apply(left: Image) = {
      MatchResult(
        left.asScala.forall((ic) => left.getValue(ic) == expectedImage.getValue(ic)),
        "Image did not look as expected.",
        "Image looked as expected"
      )
    }
  }

  def lookLike(expectedImage: Image) = new LookLikeMatcherCls(expectedImage)
}

class FindBeadsSpec extends FlatSpec with Matchers with LookLikeMatcher {

  val testCmdLineArgs = Config(segchannel = 0, segplane = 4, beadradius = 24.0)

  val imURL = new URL(
    "https://storage.googleapis.com/find-beads/test_bead_input.ome.tif")
  val imOutURL = new URL(
    "https://storage.googleapis.com/find-beads/test_bead_output.ome.tif")

  def readImageFromURL(url: URL): Image = {
    val tf = File.createTempFile("test-image", ".ome.tif")
    tf.deleteOnExit()

    val outStream = new BufferedOutputStream(new FileOutputStream(tf))
    val conn = url.openConnection()
    conn.connect()
    val inStream = new BufferedInputStream(conn.getInputStream())

    var currByte = inStream.read()

    while (currByte != -1) {
      outStream.write(currByte)
      currByte = inStream.read()
    }

    outStream.close()
    inStream.close()

    new ImageReader().read(tf.getAbsolutePath())
  }

  // TODO: download these only the first time you run the test.
  val im = readImageFromURL(imURL)
  val imOut = readImageFromURL(imOutURL)

  it should "correctly segment a bead image" taggedAs(SlowTest) in {
    FindBeads.maskFromImage(im, testCmdLineArgs) should lookLike(imOut)
  }

  it should "correctly calculate the centroids of objects in a mask" in {
    val imTest = ImageFactory.create(ImageCoordinate.createCoordXYZCT(11, 11, 1, 1, 1), 1.0f)
    val expected = Map[Float, ArrayBuffer[Double]]((1.0f -> ArrayBuffer(5.0, 5.0)))
    FindBeads.centroids(imTest) should equal(expected)

    val partialExpectedRealImage = Map[Float, ArrayBuffer[Double]](
      1.0f -> ArrayBuffer(97.11400110680687, 44.54676258992806),
      2.0f -> ArrayBuffer(49.38269550748752, 47.09317803660566))

    val realCentroids = FindBeads.centroids(imOut)

    realCentroids.get(1.0f) should equal(partialExpectedRealImage.get(1.0f))
    realCentroids.get(2.0f) should equal(partialExpectedRealImage.get(2.0f))
  }

  it should "correctly perform projection onto a vector" in {
    val origin = Array(0.0, 0.0)
    val pointOnLine = Array(3.0, 0.0)
    val pointToProject = Array(7.0, 2.0)
    val expected = Array(7.0, 0.0)

    val result = FindBeads.projectPointOntoVector(origin, pointToProject, pointOnLine)
    result should contain theSameElementsInOrderAs expected
  }

  it should "correctly calculate region borders of a Voronoi diagram" in {
    val points = Array(ArrayBuffer(-2.0, 2.0), ArrayBuffer(2.0, 2.0))
    val icOn = ImageCoordinate.createCoordXYZCT(0, 2, 0, 0, 0)
    val icOff = ImageCoordinate.createCoordXYZCT(1, 2, 0, 0, 0)
    FindBeads.isOnVoronoiBorder(points, icOff) should be(false)
    FindBeads.isOnVoronoiBorder(points, icOn) should be(true)
  }

  it should "caclulate a reasonable maximum bead size for a given radius" in {
    val rad = 10.0
    val maxSize = FindBeads.calculateMaxSizeFromRadius(rad)

    maxSize.asInstanceOf[Double] should be < (Math.pow(2*rad + 1, 2))
    maxSize.asInstanceOf[Double] should be > (Math.PI * Math.pow(rad, 2))
  }

  it should "calculate a reasonable minimum bead size for a given radius" in {
    val rad = 10.0
    val minSize = FindBeads.calculateMinSizeFromRadius(rad)

    minSize.asInstanceOf[Double] should be < (Math.PI * Math.pow(rad, 2))
    minSize.asInstanceOf[Double] should be > (0.25 * Math.PI * Math.pow(rad, 2))
  }
}







