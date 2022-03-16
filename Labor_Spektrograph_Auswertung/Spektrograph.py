
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.linalg import fractional_matrix_power

images_mit = []
images_ohne = []

_height=0
_width =0
_depth = 0

_img_nr = 0

_pp7r = 0

def showimage(im):

    plt.figure()
    plt.title("Original File Data")
    plt.imshow(im)
    plt.show()

def load_images():
    # open image
    print("\nopening image")

    intens_mit = [1/5,1/5,2,1/30,1/60,1/125,1/125,1/50,1/40,1/2,1/25,1/25,1/10,1,4]
    print("len intens mit "+ str(len(intens_mit)))
    intens_ohne = [1 / 5, 1 / 5, 1, 1, 1 / 2, 1, 1 , 1 , 2, 2, 4, 5,4, 4,4]
    print("len intens ohne " + str(len(intens_ohne)))


    for i in range(5,35):
        index_pre= "";
        if i < 10:
            index_pre = "000"
        else:
            index_pre = "00"

        im = Image.open("2020_06_15/IMG_"+index_pre+str(i)+".jpg")
        print("load 2020_06_15/IMG_"+index_pre+str(i)+".jpg")

        if(i%2 == 1):
            images_mit.append(im)
        else:
            images_ohne.append(im)

    print("len of img containers")
    print(len(images_mit))
    print(len(images_ohne))
    _img_nr = len(images_mit)

def getPixelPer7Revolutions():

    pixels1 = np.asarray(images_mit[0])
    pixels2 = np.asarray(images_mit[1])

    _height,_width,_depth = pixels1.shape

    # first image pixel position
    pixelsummax1 = 0
    pixelsum_max_index1 = 0;

    for w in range(0,_width):
        pixelsum = (pixels1[int(_height/2)][w][0])
        if(pixelsummax1<pixelsum):
            pixelsummax1 = pixelsum
            pixelsum_max_index1 = w;

    # second image pixel position
    pixelsummax2 = 0
    pixelsum_max_index2 = 0;

    for w2 in range(0, _width):
        pixelsum2 = (pixels2[int(_height / 2)][w2][0])
        if (pixelsummax2 < pixelsum2):
            pixelsummax2 = pixelsum2
            pixelsum_max_index2 = w2;

    _pp7r = pixelsum_max_index1-pixelsum_max_index2
    print("diff pp7r "+str(_pp7r))

def stickImages(list,filename):
    begin = int(500)
    end = int(begin+_pp7r)

    p3 = np.asarray(list[0]).transpose(1,0,2)
    print(p3.shape)

    for i in range(0,15):
        p1 = np.asarray(list[i]) # make img to array

        p2 = p1.transpose((1,0,2)) # transpose
        p4 = p2[begin:begin+1586] # cut out part of image
        print("stick part "+str(p4.shape)) # shape of cut out part

        p3 = np.vstack((p3, p4))    # join

    print("whole img dimensions "+str(p3.shape))
    p5 = p3.transpose(1,0,2)
    p5 = p5 * float(255 / np.max(p5))  # intensitäten angleichen
    image = Image.fromarray(np.uint8(p5))
    showimage(image)

    # save image
    image.save(filename)

    return p3;

def getSpectrum(img):
    res = 2;
    height = 200;

    matrix = np.asarray(img).transpose(1,0,2)

    intens = [] #list of intensities
    for i in range(len(matrix)):
        pixsum = 0;
        for j in range(1628, 1828): # über diese zeilen
            pixsum += np.sum(matrix[i][j])
        intens.append(pixsum)

    plt.figure()
    plt.title("Spectrum")
    plt.plot(range(len(intens)),intens)
    plt.show()




load_images() #load
getPixelPer7Revolutions() # pixel transition per 7*50um camera transition
mat1 = stickImages(images_mit,"spektrum_mit.jpg"); # get panorama
mat2 = stickImages(images_ohne,"spektrum_ohne.jpg");

getSpectrum(Image.open("spektrum_mit.jpg")); # get spectrum intensitie over pixels
getSpectrum(Image.open("spektrum_ohne.jpg"));