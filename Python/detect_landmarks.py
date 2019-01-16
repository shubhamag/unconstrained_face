#!/usr/bin/python
# Preprocess folder for face reconstruction.
# Detect faces
# Scale and crop images
# Perform landmark alignment

import csv
import dlib
import glob
import os
import scipy as sp
import scipy.io
import shutil
import skimage.io
import skimage.transform
import sys

if len(sys.argv) != 2:
    print(
        "Must provide directory with face images\n"
        "    ./detect_landmarks.py ../MattDamon\n")
    print(sys.argv)
    exit()

predictor_path = 'Python/shape_predictor_68_face_landmarks.dat'
faces_path = sys.argv[1]

detector = dlib.get_frontal_face_detector()
predictor = dlib.shape_predictor(predictor_path)

files = glob.glob(os.path.join(faces_path, "*.jpg"))
export_path = os.path.join(faces_path, 'Export')
if os.path.exists(export_path):
    shutil.rmtree(export_path)
os.mkdir(export_path)

fo = open(os.path.join(export_path, "landmarks.csv"),"wb")
writer = csv.writer(fo)

for fi, f in enumerate(files):
    img = skimage.io.imread(f)

    # Detect faces
    dets, scores, idx = detector.run(img)
    if len(dets) and scores[0] > 1.25:
        # Scale image
        d = dets[0]
        scale = 250.0 / max(d.right() - d.left(), d.bottom() - d.top())
        img = skimage.img_as_ubyte(skimage.transform.rescale(img, scale))

        # Crop image
        dets = detector(img)
        if len(dets):
            d = dets[0]
            left = max(0, d.left() - 75)
            right = min(img.shape[1], d.right() + 75)
            top = max(0, d.top() - 75)
            bottom = min(img.shape[0], d.bottom() + 75)
            img = img[top:bottom, left:right].copy()

            # Perform landmark alignment
            dets = detector(img)
            if len(dets):
                filename = os.path.join(export_path, "{:0>3d}_face.jpg".format(fi))
                skimage.io.imsave(filename, img)
                shape = predictor(img, dets[0])
                landmarks = []
                for i in range(0,68):
                    landmarks.append(shape.part(i).y)
                    landmarks.append(shape.part(i).x)
                writer.writerow(landmarks)
