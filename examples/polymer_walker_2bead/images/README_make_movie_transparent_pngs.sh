
# If you have images with alpha (transparency), you can change the background to black:
#for f in *.png; do
#  bn=`basename "$f" .png`
#  convert "$f" -background black -alpha remove "${bn}_bbk.png"
#done
#ffmpeg -r 48 -f image2 -i walker_2bead_ring24_mass1_%04d_bbk.png  -codec:v libx264 -crf 20 walker_2bead_ring24_mass1_bbk.mp4

# Otherwise, by default, transparent pixels are converted to white
ffmpeg -r 48 -f image2 -i walker_2bead_ring24_mass1_%04d.png  -codec:v libx264 -crf 20 walker_2bead_ring24_mass1_wbk.mp4
