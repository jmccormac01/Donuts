from ..donuts import Donuts


# Test the donuts code using real data
def test_full_integration():
	# initialise the class with the settings needed
	# upon generation of the reference image
	d = Donuts(refimage = 'IMAGE80520160114005507.fits', image_ext = 0, exposure = 'EXPOSURE', 
			normalise = True, subtract_bkg = True, prescan_width = 20, overscan_width = 20, boarder = 64, ntiles = 32)
	# print a summary of the setup
	d.print_summary()
	# assumes all the settings from the ref image generation
	# and calculates the shift between the images
	imlist=['IMAGE80520160114005520.fits','IMAGE80520160114005533.fits']
	for image in imlist:
		x, y = d.measure_shift(checkimage=image)


