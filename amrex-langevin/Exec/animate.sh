#!/bin/sh

# Use ffmpeg to animate plot slices
ffmpeg -r 16 -i plt00000_phi1_%05d.png -s 1440x1080 -vcodec libx264 -crf 32 -pix_fmt yuv420p phi1.mp4
ffmpeg -r 16 -i plt00000_phi2_%05d.png -s 1440x1080 -vcodec libx264 -crf 32 -pix_fmt yuv420p phi2.mp4
ffmpeg -r 16 -i plt00000_phi3_%05d.png -s 1440x1080 -vcodec libx264 -crf 32 -pix_fmt yuv420p phi3.mp4
ffmpeg -r 16 -i plt00000_phi4_%05d.png -s 1440x1080 -vcodec libx264 -crf 32 -pix_fmt yuv420p phi4.mp4

