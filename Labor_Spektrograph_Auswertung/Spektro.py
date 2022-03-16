# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 20:14:50 2020

@author: Patrik
"""

import numpy as np
from PIL import Image
import os, os.path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

def Disp_R(px):
    return (3.13e-7)*px**2 - 0.0155*px + 653.68

#Function for merging 2 pictures
def get_concat_h(im1, im2):
    w1,h1 = im1.size
    w2,h2 = im2.size
    dst = Image.new('RGB', (w1 + w2, h1))
    dst.paste(im1, (0, 0))
    dst.paste(im2, (w1, 0))
    return dst

#Function to crop image to square image from upper left corner
def Upperleftsquare(im):
    #im = imori.copy()
    w,h = im.size  # Get dimensions
    if w > h:
        return im.crop((0,0,h,h))
    if w < h:
        return im.crop((0,0,w,w))
    else:
        print('No cropping needed, Image already square')
        return

#Function to get pixelvalues as (r,g,b)-tupels from image
def GetPixelvalues(im):
    width, height = im.size
    pixel_values = list(im.getdata())
    pixel_values = np.array(pixel_values).reshape((height, width, 3))
    return pixel_values

#Function to store r,g,b values in seperate matrices from pixelvalues((r,g,b)-tupels)
def Getrgb(pixel_values):
    h = len(pixel_values)
    b = len(pixel_values[0])
    Red = np.zeros((h,b))
    Green = np.copy(Red)
    Blue = np.copy(Red)
    for i in range(h):
        for j in range(b):
            Red[i][j] = np.copy(pixel_values[i][j][0])
            Green[i][j] = np.copy(pixel_values[i][j][1])
            Blue[i][j] = np.copy(pixel_values[i][j][2])
    return Red, Green, Blue

#Funktion im Intensitaeten von den Pixel-Werten zu erhalten
def GetIntensities(pixel_values):
    h = len(pixel_values)
    b = len(pixel_values[0])
    Intensities = np.zeros((h, b))
    for i in range(h):
        for j in range(b):
            Intensities[i][j] = pixel_values[i][j][0]+pixel_values[i][j][1]+\
            pixel_values[i][j][2]
    return Intensities

#Umdrehungen am Rad beim Aufnehmen des SPektrums
steps = [7,7,7,7,7,7,7,7,7,7,7,7,7]
imgs_ref = []
imgs_hal = []

#print 'steps',len(steps)

#Importieren der Bilder
path_hal = "2020_06_15/absorb"
valid_images = [".jpg",".gif",".png",".tga"]
for f in os.listdir(path_hal):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_images:
        continue
    imgs_hal.append(Image.open(os.path.join(path_hal,f)))

path_ref = "2020_06_15/ref"
valid_images = [".jpg",".gif",".png",".tga"]
for f in os.listdir(path_ref):
    ext = os.path.splitext(f)[1]
    if ext.lower() not in valid_images:
        continue
    imgs_ref.append(Image.open(os.path.join(path_ref,f)))
#
# for i in range(len(imgs_ref)):
#     imgplot = plt.imshow(imgs_ref[i])
#     plt.show()
#
# for i in range(len(imgs_hal)):
#     imgplot = plt.imshow(imgs_hal[i])
#     plt.show()

#Wegschneiden des Teils der Bilder wo nichts aufgenommen wird
imgs_ref_crp = []
imgs_hal_crp = []
#Cropping images
for i in range(len(imgs_ref)):
    imgs_ref_crp.append(Upperleftsquare(imgs_ref[i]))
for i in range(len(imgs_ref)):
    imgs_hal_crp.append(Upperleftsquare(imgs_hal[i]))

# pix_ref_0 = GetPixelvalues(imgs_ref[0])
# rgb_ref_0 = Getrgb(pix_ref_0)
# pix_ref_1 = GetPixelvalues(imgs_ref[1])
# rgb_ref_1 = Getrgb(pix_ref_1)
#
# max_ref_0 = np.argmax(rgb_ref_0[0][0])
# max_ref_1 = np.argmax(rgb_ref_1[0][0])

#print max_ref_0,max_ref_1

#delta_pix_p_step = abs(max_ref_1-max_ref_0)/steps[0]
delta_pix_p_step = 219
#print delta_pix_p_step

img_sprectrum_ref = imgs_ref_crp[0].copy()

#Zusammenfuegen Referenz
for i in range(len(steps)):
    w, h = imgs_ref_crp[i+1].size
    #print w-delta_pix_p_step*steps[i],0,w,h
    croped_im = imgs_ref_crp[i+1].crop((w-delta_pix_p_step*steps[i]-100,0,w,h))
    #imgplot = plt.imshow(croped_im)
    # plt.show()
    img_sprectrum_ref = get_concat_h(img_sprectrum_ref,croped_im)

#Zusammenfuegen Halogen
img_sprectrum_hal = imgs_hal_crp[0].copy()
for i in range(len(steps)):
    w, h = imgs_hal_crp[i+1].size
    #print w-delta_pix_p_step*steps[i],0,w,h
    croped_im = imgs_hal_crp[i+1].crop((w-delta_pix_p_step*steps[i]-100,0,w,h))
    #imgplot = plt.imshow(croped_im)
    #plt.show()
    img_sprectrum_hal = get_concat_h(img_sprectrum_hal,croped_im)

#img_sprectrum_ref.save('Spektrum_Referenz.jpg')
#img_sprectrum_hal.save('Spektrum_Halogen.jpg')
# imgplot = plt.imshow(img_sprectrum_ref)
# plt.show()

#Berechenung der Intensitaeten Referenz Spektrum
w_ref_spec, h_ref_spec = img_sprectrum_ref.size
spectrum_ref_rgb = GetIntensities(GetPixelvalues(img_sprectrum_ref.crop\
                    ((0,int(h_ref_spec/2)-200,w_ref_spec,int(h_ref_spec/2)+200))))

inten_vector_ref = np.sum(spectrum_ref_rgb,axis=0)

xvalues = np.arange(0,len(inten_vector_ref))
plt.figure(dpi=800)
plt.plot(xvalues,inten_vector_ref, linewidth = 0.7)
plt.xlabel('Pixelposition')
plt.ylabel('I / a.u.')
plt.title('Referenz')
plt.savefig('Ref.png')
plt.show()

#Berechnung der Intensitaeten Halogen Spektrum
w_hal_spec, h_hal_spec = img_sprectrum_hal.size
spectrum_hal_rgb = GetIntensities(GetPixelvalues(img_sprectrum_hal.crop\
                    ((0,int(h_hal_spec/2)-200,w_hal_spec,int(h_hal_spec/2)+200))))

inten_vector_hal = np.sum(spectrum_hal_rgb,axis=0)

xvalues = np.arange(0,len(inten_vector_hal))
plt.figure(dpi=800)
plt.plot(xvalues,inten_vector_hal, linewidth = 0.7)
plt.xlabel('Pixelposition')
plt.ylabel('I / a.u.')
plt.title('Absorption')
plt.savefig('Hal.png')
plt.show()

#Maxima Position fur Berechnung Dispersionsrel
print( 'Rot',np.argmax(inten_vector_ref[:5000]))
print( 'Gelb',5000+np.argmax(inten_vector_ref[5000:9000]),5000+np.argmax(inten_vector_ref[5000:4950+np.argmax(inten_vector_ref[5000:9000])]))
print( 'Hellgruen',9000+np.argmax(inten_vector_ref[9000:14000]))
print( 'Gruen',14000+np.argmax(inten_vector_ref[14000:18000]))
print( 'Blaugruen',16000+np.argmax(inten_vector_ref[16000:18000]))
print( 'Blau',19000+np.argmax(inten_vector_ref[19000:20000]))
print( 'Violett',20000+np.argmax(inten_vector_ref[20000:25000]))
#print( 'Dunkelviolett',25000+np.argmax(inten_vector_ref[25000:35000]))

xvalues = Disp_R(np.arange(0, len(inten_vector_ref)))
plt.figure(dpi=800)
plt.plot(xvalues,inten_vector_ref, linewidth = 0.7)
plt.xlim(max(xvalues),min(xvalues))
plt.xlabel('WL / nm')
plt.ylabel('I / a.u.')
plt.title('Wellenlängen Referenz')
plt.savefig('Ref_lambda.png')
plt.show()

xvalues = Disp_R(np.arange(0, len(inten_vector_hal)))
plt.figure(dpi=800)
plt.plot(xvalues, inten_vector_hal, linewidth = 0.7)
plt.xlim(max(xvalues), min(xvalues))
plt.xlabel('WL / nm')
plt.ylabel('I / a.u.')
plt.title('Wellenlängen Absorption')
plt.savefig('Hal_lambda.png')
plt.show()

#Halbwertsbreite des Gruen Peaks
# xvalues = Disp_R(np.arange(0,len(inten_vector_ref)))
# plt.plot(xvalues,inten_vector_ref)
# plt.xlim(549.5,550)
# plt.xlabel('WL / nm')
# plt.ylabel('I / a.u.')
# plt.grid()
# plt.savefig('Ref_lambda.png')
# plt.show()

#Verharltnis wellelaenge zu pixel
print(max(xvalues)/len(inten_vector_ref))

#Wellenlaengen der Bandkanten
# xvalues = np.arange(0,len(inten_vector_hal))
# plt.figure(dpi=1200)
# plt.plot(xvalues,inten_vector_hal)
# plt.xlim(12000,14000)
# plt.xlabel('WL / nm')
# plt.ylabel('I / a.u.')
# plt.grid()
# plt.savefig('lambda_bandkant.png')
# plt.show()

#Positionen der Minima fue Dissozoationsenergie
minimums = []
for i in range(10):
    if i == 5:
        min = 9510 + i*200+50 + np.argmin(inten_vector_hal[9510 + i*200 + 50:9510 + (i+1)*200 + 50])
    else:
        min = 9510 + i * 200 + np.argmin(inten_vector_hal[9510 + i * 200:9510 + (i + 1) * 200])
    #print min
    minimums.append(Disp_R(min))

print(minimums)