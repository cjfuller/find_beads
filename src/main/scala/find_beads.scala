package io.cjf.findbeads

import java.io.{File, PrintWriter}
import scala.collection.JavaConverters._
import scala.collection.mutable.ArrayBuffer
import util.control.Breaks._

import edu.stanford.cfuller.imageanalysistools.filter.{
  Filter,
  LabelFilter,
  MaskFilter,
  MaximumSeparabilityThresholdingFilter,
  SizeAbsoluteFilter
}
import edu.stanford.cfuller.imageanalysistools.frontend.{
  LocalAnalysis,
  MaximumIntensityProjection
}
import edu.stanford.cfuller.imageanalysistools.image.{
  Histogram,
  Image,
  ImageCoordinate,
  ImageFactory,
  ImageSet,
  WritableImage
}
import edu.stanford.cfuller.imageanalysistools.image.io.{
  ImageReader,
  ImageWriter
}
import edu.stanford.cfuller.imageanalysistools.meta.parameters.{
  Parameter,
  ParameterDictionary
}
import edu.stanford.cfuller.imageanalysistools.metric.IntensityPerPixelMetric

import io.cjf.findbeads.config.Config


object FindBeads {
  implicit def javaToScalaDouble(d: java.lang.Double) = d.doubleValue

  /**
    * Calculate the minimum allowed size of a bead from the supplied radius.
    *
    * This is set to be slightly smaller than a third of a circle of that
    * radius. Note that this is smaller than any of the returned regions should
    * be, but making the cutoff this small is useful for dividing up clumps of
    * beads where several rounds of recursive thresholding may make the regions
    * quite small temporarily.
    *
    * @param [Double] radius: the radius of the bead in pixels.
    */
  def calculateMinSizeFromRadius(radius: Double): Int = {
    Math.round(0.97 * Math.pow(radius + 1, 2)).asInstanceOf[Int]
  }

  /**
    * Calculate the maximum allowed size of a bead from the supplied radius.
    *
    * This is set to be slightly larger than a circle of that radius.
    *
    * @see #calculateMinSizeFromRadius
    */
  def calculateMaxSizeFromRadius(radius: Double): Int = {
    Math.round(3.2 * Math.pow(radius + 1, 2)).asInstanceOf[Int]
  }

  /**
    * Thresholds the values in an array according to Otsu's method (Otsu,
    * 1979, DOI: 10.1109/TSMC.1979.4310076).
    *
    * @param [Array[Int]] values: an array (1-D) containing the values to be thresholded.
    * @return [Int] the threshold value; values <= to this number are below the threshold.
    */
  def graythresh(values: Array[Float]): Float = {
    val numSteps = 1000.0
    val min = values.min.asInstanceOf[Double]
    val max = values.max.asInstanceOf[Double]
    if (max == min) {
      return min.asInstanceOf[Float]
    }
    val mean: Double = values.sum / values.length
    val sortedValues = values.sorted
    val stepSize = (max - min) / numSteps

    var bestEta = 0.0
    var bestThresholdValue = 0.0
    var countsUptoK = 0
    var mu = 0.0

    (min to max by stepSize).foreach((k: Double) => {
      val rest = sortedValues.drop(countsUptoK)
      val inSegment = rest.takeWhile((v: Float) => v <= k)
      countsUptoK += inSegment.length
      mu += inSegment.sum / sortedValues.length
      val omega: Double = countsUptoK.asInstanceOf[Double] / sortedValues.length

      if (omega > 0 && omega < 1) {
        val eta = omega * (1 - omega) * Math.pow((mean - mu)/(1 - omega) - mu/omega, 2)
        if (eta > bestEta) {
            bestEta = eta
            bestThresholdValue = k
        }
      }
    })
    bestThresholdValue.asInstanceOf[Float]
  }

  /**
  * Recursively thresholds regions in a supplied mask and image using the
  * method described in Xiong et al. (DOI: 10.1109/ICIP.2006.312365).
  *
  * @param [ParameterDictionary] params: a ParameterDictionary specifying
  *   max_size and min_size parameters, which control the maximum size of
  *   regions before they are recursively thresholded to break them up, and the
  *   minimum size of regions before they are discarded.
  * @param [Image] im: the original image being segmented. This will not be
  *   modified.
  * @param [Image] mask: the initial segmentation mask for the supplied image.
  *   Regions in this mask may be divided up or discarded.
  */
  def recursiveThreshold(params: ParameterDictionary, im: Image, mask: WritableImage) {
    val h = new Histogram(mask)
    var changed = false
    var discardList: Set[Int] = Set()

    val minSize = Integer.parseInt(params.getValueForKey("min_size"))
    val maxSize = Integer.parseInt(params.getValueForKey("max_size"))

    val histogramRange = 1 to h.getMaxValue()
    val coordinatesByMaskValue: scala.collection.mutable.Map[Int, List[ImageCoordinate]] = (
      scala.collection.mutable.Map())
    var thresholdsByMaskValue: Map[Int, Float] = Map()

    im.asScala.foreach {(ic) => {
      val maskVal = mask.getValue(ic).asInstanceOf[Int]
      if (h.getCounts(maskVal) > maxSize && maskVal > 0) {
        val currList = coordinatesByMaskValue.getOrElse(maskVal, List[ImageCoordinate]())
        coordinatesByMaskValue(maskVal) = ImageCoordinate.cloneCoord(ic) :: currList
      }
    }}

    coordinatesByMaskValue.foreach {case (maskVal, imCoords) => {
      val imValues = imCoords.map((ic: ImageCoordinate) => im.getValue(ic))
      val threshold = graythresh(imValues.toArray)
      thresholdsByMaskValue = thresholdsByMaskValue + (maskVal -> threshold)
    }}

    histogramRange.foreach((i: Int) => {
      if (h.getCounts(i) > 0 && h.getCounts(i) < minSize) {
        discardList = discardList + i
      }
    })

    im.asScala.foreach((ic) => {
      val maskVal = mask.getValue(ic).asInstanceOf[Int]
      thresholdsByMaskValue.get(maskVal) match {
        case Some(thresh) => {
          if (im.getValue(ic) <= thresh) {
            mask.setValue(ic, 0)
            changed = true
          }
        }
        case None => (
          if (discardList contains mask.getValue(ic).asInstanceOf[Int]) {
              mask.setValue(ic, 0)
              changed = true
          }
        )
      }
    })

    if (changed) {
      val lf = new LabelFilter()
      lf.apply(mask)
      recursiveThreshold(params, im, mask)
    }
  }

  /**
    * Generate a segmented mask of beads froom an image.
    *
    * @param [Image] im: the image to segment
    * @param [Config] opts: the command-line options
    */
  def maskFromImage(im: Image, opts: Config): Image =  {
    val radius = opts.beadradius
    val minSize = calculateMinSizeFromRadius(radius)
    val maxSize = calculateMaxSizeFromRadius(radius)

    val sizes = ImageCoordinate.cloneCoord(im.getDimensionSizes())
    sizes.set(ImageCoordinate.C, 1)
    sizes.set(ImageCoordinate.Z, 1)

    val im0 = ImageCoordinate.createCoordXYZCT(0, 0, 0, 0, 0)
    im0.set(ImageCoordinate.C, opts.segchannel)
    im0.set(ImageCoordinate.Z, opts.segplane)
    val imToSegment = im.subImage(sizes, im0).getWritableInstance()

    val params = ParameterDictionary.emptyDictionary()
    params.setValueForKey("min_size", minSize.toString())
    params.setValueForKey("max_size", maxSize.toString())

    val imCopy = ImageFactory.createWritable(imToSegment)

    val lf = new LabelFilter()

    val filters: Array[Filter] = Array(
      new MaximumSeparabilityThresholdingFilter(),
      lf
    )

    filters.foreach((f) => {
      f.setParameters(params)
      f.setReferenceImage(imCopy)
      f.apply(imToSegment)
    })
    recursiveThreshold(params, imCopy, imToSegment)

    val saf = new SizeAbsoluteFilter()
    saf.setParameters(params)
    saf.apply(imToSegment)
    val cens = centroids(imToSegment)

    val finalMask = ImageFactory.createWritable(imToSegment)

    imToSegment.asScala.foreach((ic: ImageCoordinate) => {
      finalMask.setValue(ic, 0.0f)
      cens.foreach { case (region, cen) => {
        if (Math.hypot(cen(0) - ic.get(0), cen(1) - ic.get(1)) <= radius) {
          finalMask.setValue(ic, region)
        }
      }}
    })

    finalMask.asScala.foreach((ic: ImageCoordinate) => {
      if (finalMask.getValue(ic) > 0 && isOnVoronoiBorder(cens.values.toSeq, ic)) {
        finalMask.setValue(ic, 0)
      }
    })

    lf.apply(finalMask)
    saf.apply(finalMask)
    lf.apply(finalMask)

    finalMask
  }

  /**
    * Vector L2 norm
    */
  def vnorm(v: Seq[Double]): Seq[Double] = {
    val len = Math.sqrt(v.map((vi) => Math.pow(vi, 2)).sum)
    v.map((vi) => vi / len)
  }

  /**
    * Vector subtraction
    */
  def vsub(v1: Seq[Double], v0: Seq[Double]): Seq[Double] = {
    v1.zip(v0).map{case (vi, vj) => vi - vj}
  }

  /**
    * Vector addition
    */
  def vadd(v1: Seq[Double], v0: Seq[Double]): Seq[Double] = {
    v1.zip(v0).map{case (vi, vj) => vi + vj}
  }

  /**
    * Vector inner product
    */
  def vdot(v0: Seq[Double], v1: Seq[Double]): Double = {
    v1.zip(v0).map{case (vi, vj) => vi * vj}.sum
  }

  /**
    *  Scalar on vector multiplication
    */
  def stimes(s: Double, v: Seq[Double]): Seq[Double] = {
    v.map((vi) => s * vi)
  }

  /**
    * Projects a point in space onto a specified line
    *
    * @param [Seq] origin the point in space the will serve as the origin
    *  for the purposes of the projection (this should be on the line)
    * @param [Seq] point_to_project the point in space that will be projected
    *  this should be in the same coordinate system in which the origin is
    *  specified, not relative to the origin
    * @param [Seq] point_on_line another point on the line specified in the
    *  same coordinate system in which the origin is specified, not relative to
    *  the origin
    *
    * @return [Seq] the projected point (in the same coordinate system as
    *  the other points were specified)
    */
  def projectPointOntoVector(
    origin: Seq[Double], pointToProject: Seq[Double],
    pointOnLine: Seq[Double]): Seq[Double] = {

    val unitVec = vnorm(vsub(pointOnLine, origin))
    val projVec = vsub(pointToProject, origin)
    vadd(stimes(vdot(projVec, unitVec), unitVec), origin)
  }

  /**
    * Checks if a given coordinate would be approximately on the boundary
    * between two regions of a Voronoi diagram of constructed from a set of
    * points.  The approximation is calculated such that no two regions in the
    * Voronoi diagram would be 8-connected.
    *
    * @param [Seq] points: a seq of two-element arrays containing the x,y
    * coordinates of the points on which the diagram is calculated.
    * @param [ImageCoordinate] ic: an ImageCoordinate specifying the location to check.
    * @return [Boolean] whether the point is on the border between two regions.
    */
  def isOnVoronoiBorder(points: Seq[ArrayBuffer[Double]], ic: ImageCoordinate): Boolean = {
    val x = ic.get(ImageCoordinate.X)
    val y = ic.get(ImageCoordinate.Y)

    var closestIndex = 0
    var nextIndex = 0
    var closestDist = Double.MaxValue
    var nextDist = Double.MaxValue

    points.zip(points.indices).foreach { case (p, i) => {
      val dist = Math.hypot(p(0) - x, p(1) - y)

      if (dist < closestDist) {
        nextDist = closestDist
        nextIndex = closestIndex
        closestDist = dist
        closestIndex = i
      } else if (dist < nextDist) {
        nextDist = dist
        nextIndex = i
      }
    }}

    val projPoint = projectPointOntoVector(
      points(closestIndex), Array(x, y), points(nextIndex))

    val nextDistProj = Math.hypot(
      points(nextIndex)(0) - projPoint(0),
      points(nextIndex)(1) - projPoint(1))

    val closestDistProj = Math.hypot(
      points(closestIndex)(0) - projPoint(0),
      points(closestIndex)(1) - projPoint(1))

    val cutoff = 1.01 * Math.sqrt(2)

    (nextDistProj - closestDistProj < cutoff)
  }

  def centroids(mask: Image): Map[Float, ArrayBuffer[Double]] = {
    var cens = Map[Float, ArrayBuffer[Double]]()

    mask.asScala.foreach((ic: ImageCoordinate) => {
      if (mask.getValue(ic) > 0) {
        val maskValue = mask.getValue(ic)
        val cen = cens.getOrElse(maskValue, ArrayBuffer[Double](0.0, 0.0))
        val x: Int = ic.get(0)
        val y: Int = ic.get(1)
        cen(0) += x
        cen(1) += y
        cens = cens + ((maskValue, cen))
      }
    })

    val h = new Histogram(mask)

    cens.foreach { case (k, cen) => {
      val count = h.getCounts(k.asInstanceOf[Int])
      cen(0) /= count
      cen(1) /= count
    }}
    cens
  }

  /**
    * Write the output data and mask to files.
    *
    * @param fnOrig: the filename of the image being processed
    * @param quantStr: the quantification data to write
    * @param mask: the mask to write to file
    */
  def writeOutput(fnOrig: File, quantStr: String, mask: Image) {
    val dir = fnOrig.getParent()
    val base = fnOrig.getName().replaceAll(".ome.tif", "")

    val maskDir = new File(dir, "output_mask")
    val quantDir = new File(dir, "quantification")
    maskDir.mkdir()
    quantDir.mkdir()

    val maskExt = "_mask.ome.tif"
    val quantExt = "_quant.txt"

    mask.writeToFile((new File(maskDir, base + maskExt)).getAbsolutePath())
    val pw = new PrintWriter(new File(quantDir, base + quantExt))
    pw.println(quantStr)
    pw.close()
  }

  /**
    * Process a single file, which consists of creating a mask, quantifying
    * regions, and writing output.
    *
    * @param [File] f: the image file to process
    * @param [Config] config: the configuration object of command line options.
    */
  def processFile(f: File, config: Config) {
    println(s"Processing ${f.getName()}")
    val ir = new ImageReader()
    val im = ir.read(f.getAbsolutePath());
    val mask = FindBeads.maskFromImage(im, config)
    val proj = MaximumIntensityProjection.projectImage(im)
    val ims = proj.splitChannels()
    val is = new ImageSet(ParameterDictionary.emptyDictionary())

    ims.asScala.foreach {(imc) => is.addImageWithImage(imc)}

    val metric = new IntensityPerPixelMetric()
    val quant = metric.quantify(mask, is)
    val outdat = LocalAnalysis.generateDataOutputString(quant, null)
      .trim().replaceAll(" ", ",")
    writeOutput(f, outdat, mask)
  }

}
