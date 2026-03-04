#! python

'''
This demo extracts all images of a PDF as PNG files, whether they are
referenced by pages or not.
It scans through all objects and selects /Type/XObject with /Subtype/Image.
So runtime is determined by number of objects and image volume.
Usage:
extract_img2.py input.pdf
'''

from __future__ import print_function
import fitz
import sys, time, re

checkXO = r"/Type(?= */XObject)"  # finds "/Type/XObject"
checkIM = r"/Subtype(?= */Image)"  # finds "/Subtype/Image"

t0 = time.process_time()
doc = fitz.open("PDFs/Sushi.pdf")
imgcount = 0
lenXREF = doc._getXrefLength()  # number of objects - do not use entry 0!

# display some file info
print("file: %s, pages: %s, objects: %s" % (sys.argv[1], len(doc), lenXREF - 1))

for i in range(1, lenXREF):  # scan through all objects
    text = doc._getXrefString(i)  # string defining the object
    isXObject = re.search(checkXO, text)  # tests for XObject
    isImage = re.search(checkIM, text)  # tests for Image
    # if not isXObject or not isImage:  # not an image object if not both True
    #     continue
    imgcount += 1
    try:
        pix = fitz.Pixmap(doc, i)  # make pixmap from image
    except RuntimeError:  # New code added!
        pix = None
        continue
    if pix.n < 5:  # can be saved as PNG
        pix.writePNG("img-%s.png" % (i,))
    else:  # must convert the CMYK first
        pix0 = fitz.Pixmap(fitz.csRGB, pix)
        pix0.writePNG("img-%s.png" % (i,))
        pix0 = None  # free Pixmap resources
    pix = None  # free Pixmap resources

t1 = time.process_time()
print("run time", round(t1 - t0, 2))
print("extracted images", imgcount)

for i in range(1, lenXREF):  # scan through all objects
    text = doc._getXrefString(i)
    isXObject = re.search(checkXO, text)  # tests for XObject
    isImage = re.search(checkIM, text)  # tests for Image
    if isXObject or isImage:
        print("##########Pic#############")
        print(i)
        #print(text)
    imgcount += 1
    try:
        pix = fitz.Pixmap(doc, i)  # make pixmap from image
    except RuntimeError:  # New code added!
        pix = None
        continue
    if pix.n < 5:  # can be saved as PNG
        pix.writePNG("img-%s.png" % (i,))
    else:  # must convert the CMYK first
        pix0 = fitz.Pixmap(fitz.csRGB, pix)
        pix0.writePNG("img-%s.png" % (i,))
        pix0 = None  # free Pixmap resources
    pix = None  # free Pixmap resources
