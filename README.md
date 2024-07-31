# Color-Based Object Detection

This project focuses on using color representations and histograms to detect and localize objects in images. The goal is to understand color spaces and data structures commonly used in color analysis. The program detects objects in an image based on color matching with a given library of object images.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Implementation Details](#implementation-details)
  - [Creating Histograms](#creating-histograms)
  - [Histogram Matching](#histogram-matching)
  - [Object Localization](#object-localization)
- [Working](#working)

## Overview
This program reads an input image containing multiple colorful objects and a library of object images. It detects and localizes objects from the library within the input image using color histograms.

## Features
- Color histogram creation for object images and input images.
- Histogram matching for object detection.
- Object localization by highlighting detected objects or drawing bounding boxes around them.
- No external libraries or scripting environments used.

## Requirements
- A Java compiler.
- A system capable of displaying images

## Installation
- Clone the repository: bgit clone https://github.com/yourusername/color-object-detection.git
- cd color-object-detection

## Usage
YourProgram.exe InputImage.rgb object1.rgb object2.rgb ... objectn.rgb
InputImage.rgb: The input image file containing multiple objects.
object1.rgb, object2.rgb, ..., objectn.rgb: The list of object images in the library.

## Examples
javac ImageDisplay.java; java ImageDisplay "testImages\multi_object_test_new\multi_object_test_new\update_rgb\Pikachu_and_Oswald_v2.rgb" "testImages\dataset\dataset\data_sample_rgb\Oswald_object.rgb" "testImages\dataset\dataset\data_sample_rgb\pikachu_object.rgb"

## Implementation Details
#### Creating Histograms
- Object Image Histograms: For each object image, create histograms for each color channel (R, G, B). Use the green chroma background to isolate object pixels from the background.
- Input Image Histogram: Create a histogram for the input image, considering it contains objects. Convert the image to HSV and create a histogram for the H (hue) channel, ignoring S (saturation) and V (value).

#### Histogram Matching
- Compare the hue histogram of the input image with the hue histograms of the object images.
- Implement a matching algorithm to detect the presence of object histograms in the input image histogram.

#### Object Localization
- Once an object is detected, highlight the pixels in the input image that match the object.
- Display a bounding box around the detected object to indicate its location.

## Working
![image](https://github.com/user-attachments/assets/b982d002-cb69-4489-bf9b-164667be76f2)
