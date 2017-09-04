
from PIL import Image
import numpy as np
import cv2
from matplotlib import pyplot as plt
def find_element_in_list(element, list_element):
    try:
        index_element = list_element.index(element)
        return index_element
    except ValueError:
        return None
bodydetection = cv2.CascadeClassifier('C:/Python27/opencv/build/etc/haarcascades/haarcascade_fullbody.xml')


face_cascade = cv2.CascadeClassifier('C:/Python27/opencv/build/etc/haarcascades/haarcascade_frontalface_default.xml')
eye_cascade = cv2.CascadeClassifier('C:/Python27/opencv/build/etc/haarcascades/haarcascade_eye.xml')
img = cv2.imread('smallface3.jpg')

org = cv2.imread('seed.jpg')
height, width = img.shape[:2]
height1, width1 = org.shape[:2]
img  = cv2.resize(img,(2*width, 2*height), interpolation = cv2.INTER_CUBIC)
org  = cv2.resize(org,(2*width1, 2*height1), interpolation = cv2.INTER_CUBIC)
detectimg=img
im2=img
im3=img

gray = cv2.cvtColor(detectimg, cv2.COLOR_BGR2GRAY)
faces = face_cascade.detectMultiScale3(gray,
    1.02,
    2,
    400,
    maxSize=(500,500),
    outputRejectLevels=True)

rects = faces[0]
neighbours = faces[1]


for (x,y,w,h) in faces[0]:
    roi_gray = gray[y:y+h, x:x+w]
    roi_color = img[y:y+h, x:x+w]
    eyes = eye_cascade.detectMultiScale(roi_gray)
    for (ex,ey,ew,eh) in eyes:
     cv2.rectangle(roi_color,(ex,ey),(ex+ew,ey+eh),(0,255,0),2)
cv2.imshow('img',img)

mask = np.zeros(img.shape[:2],np.uint8)
bgdModel = np.zeros((1,65),np.float64)
fgdModel = np.zeros((1,65),np.float64)
rect = (399,277,102,102)
cv2.grabCut(img,mask,rect,bgdModel,fgdModel,5,cv2.GC_INIT_WITH_RECT)
mask2 = np.where(((mask==0  ) | (mask==2)),0,1).astype('uint8')
img = img*mask2[:,:,np.newaxis]




body = bodydetection.detectMultiScale(gray, 1.1, 5 ,60)

grown=0
prev = 0
curr = 0 
grow=0
grow1=grow
flag=1
for y in faces[0] :
    for x in faces[2] : 
        temp = y
        prev = x
        curr = x

        grow=0
        grow1=grow
        while  prev <= curr :
            print('loop begining :',curr)
            print('loop begining :',prev , grow)
            rect = temp
            flag=0
            grow=grow+10
            grow1=grow1+10
            rect[0]=rect[0]
            rect[1]=rect[1]
            rect[2]=rect[2]+grow1
            rect[3]=rect[3]+grow
            masking = np.zeros(org.shape[:2],np.uint8)
            bgdModel1 = np.zeros((1,65),np.float64)
            fgdModel1 = np.zeros((1,65),np.float64)
            rect=(rect[0],rect[1],rect[2],rect[3])
            print(faces)
            cv2.grabCut(org,masking,rect,bgdModel1,fgdModel1,5,cv2.GC_INIT_WITH_RECT)
            maskgrow = np.where(((masking==0  )),0,1).astype('uint8')
            modimg = org*maskgrow[:,:,np.newaxis]
            grown=cv2.add(cv2.bitwise_xor(modimg,detectimg,mask=maskgrow),detectimg)
            cv2.imshow('img',grown)
            height = maskgrow.shape
            print height
            gray1 = cv2.cvtColor(grown, cv2.COLOR_BGR2GRAY)
            facegrow = face_cascade.detectMultiScale3(gray1,1.02,2,400,maxSize=(500,500),outputRejectLevels=True)
            prev=curr
            i=0
            j=0
            wlist=[]
            pixlist=[]
            for x,y in zip(facegrow[0],facegrow[2]) :
                  i=i+1
                  temp=x
                  curr=y
                  wlist.append(y)
                  pixlist.append(x)
            print(curr)
            print(prev)
            
           
                        


cv2.waitKey(0)
cv2.destroyAllWindows()



