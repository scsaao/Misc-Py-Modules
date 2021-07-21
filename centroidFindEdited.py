'''
Pratheek Nagaraj
July 1, 2011
Image Processing and Centriod Project
2. Programming Actvity Part (c)
This program will process an image file and find the centroid of the object
'''

#Import Packages
#from visual import *
from numpy import *
from math import *
import pyfits
import numdisplay


def centroidFunc(IMAGENAME,xPos,yPos,radius):
    global objectImage
    objectImage = pyfits.getdata(IMAGENAME)

    #Display the Object Image
    'Open the object image and display using DS9'
    #numdisplay.display(objectImage)

    #normalizedObject = modify( flatImage )
    normalizedObject=objectImage
    #image = noise( normalizedObject )
    image=objectImage
    #Bright Pixel
    'Get bright pixel location'
    #xPos = input("Please enter the x value for the location of the bright pixel: ")
    #yPos = input("Please enter the y value for the location of the bright pixel: ")
    print "The radius corresponds to the numer of points to be used in the calculation i.e. 1->9 and 3->49"
    #radius = input("Please enter the radius for the location of the bright pixel: ")

    brightLoc = image[ yPos-radius-1:yPos+radius, xPos-radius-1:xPos+radius ]
    print brightLoc	
    #Find Centroid
    'Call Centroid Function to locate the given centroid'
    centroidArray = calcCentroid( brightLoc )
    centroid = ( xPos + centroidArray[0], yPos + centroidArray[1] )
    uncertainty = centroidArray[2]

    #Output
    print "The centroid of the object is: " + str(centroid)
    print "The uncertainty of the calculation is: " + str(uncertainty)
    return xPos + centroidArray[0],yPos + centroidArray[1], centroidArray[2]
    'Give choice for a visual'
    #choice = raw_input("Would you like to see a visual (Y/N)? ")
    #if choice == "y" or choice == "Y" or choice == "Yes" or choice == "yes":
    #    visualFunc( brightLoc, centroidArray[0], centroidArray[1], uncertainty )
    #else:
    #    print "Goodbye!"
    
def calcCentroid( array ):
    #Total Sum
    'Use Python command to sum all elements in the array'
    sumArray = sum(array);
    
    #X Coordinate
    'Create the sum of the X Weighted Values'
    xSum = 0.0

    'Loop through the array and weight the X Values'
    weights = range( -len(array)/2+1, len(array)/2+1 )
    for element in array:
        pos = 0
        for point in element:
            xSum = xSum + weights[pos] * point
            pos = pos + 1

    'Find the X Coordinate by dividing by the total sum'
    
    xCoor = xSum / sumArray
    xCoor = round(xCoor,3)

    print xCoor
    print sumArray
    #Y Coordinate
    'Create the sum of the Y Weighted Values'
    ySum = 0.0

    'Loop through the array and weight the Y Values'
    weight = len(array)/2
    for element in array:
        ySum = ySum + weight*sum(element)
        weight = weight - 1

    'Find the Y Coordinate by dividing by the total sum'
    yCoor = ySum / sumArray
    yCoor = round(yCoor,3)

    uncertainty = 0.0
    for element in array:
        for element2 in element:
            uncertainty = uncertainty + 2 * sqrt(element2)

    uncertainty = uncertainty / sumArray
    uncertainty = round( uncertainty, 3 )

    return xCoor, yCoor, uncertainty

def modify( flatImage ):
    #Flat Image
    'Normalize the flat image using the mean'
    mean = flatImage.mean()
    normalizedFlat = flatImage / mean

    #Normalize Object Image
    'Normalize the object image using the normalized flat'
    normalizedObject = objectImage / normalizedFlat
    return normalizedObject

def noise( normalizedObject ):
    #Noise Floor
    'Get boundaties and adjust image'
    x1Bound = input("Please enter left x bound for the noise floor: ")
    x2Bound = input("Please enter right x bound for the noise floor: ")
    y1Bound = input("Please enter lower y bound for the noise floor: ")
    y2Bound = input("Please enter upper y bound for the noise floor: ")
    print "---------------------------------------------"

    subImage = normalizedObject[ y1Bound - 1: y2Bound, x1Bound - 1: x2Bound ]
    meanValue = subImage.mean()

    'New image modified with floor'
    image = normalizedObject - meanValue
    image[where(image < 0)]=0
    return image

#def visualFunc( array, x, y, unc ):
    #Display a visual of the centroid calculation
    #    'Create brightness circles with nested for loop'

    #maxValue = array.max()
    #   print maxValue
    
    #   posY = -len(array)/2
    #for row in array:
    #    posX = -len(array)/2
    #    for element in row:
    #        sphere( radius = element/(maxValue*2), pos = (posX, posY), color = color.blue )
    #        posX = posX + 1
    #    posY = posY + 1

#'Create circle for centroid'
#   sphere(radius = unc, pos = ( x, y, 1 ), color = color.red )
    

#Input
#'Open the object image and the flat image'
#objectImage = pyfits.getdata('Addedset1Averaged.fits')
#global flatImage
#flatImage = pyfits.getdata('Addedset1Averaged.fits')
#IMAGENAME="AS1P01_010T01_9000000318sxtPC00_level2.img"
#centroidFunc(IMAGENAME,298.0,289.,30.)
