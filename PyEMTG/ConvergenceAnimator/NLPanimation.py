#generator for NLP movies
def make_movie(movie_type):
    import logging
    import os
    import matplotlib
    from matplotlib.font_manager import findfont, FontProperties
        
    # Make sure imageio and PIL are actually installed
    try:
        import imageio
        from PIL import Image, ImageFont, ImageDraw
    except:
        # They arent, cant do this
        logging.info("imageio and/or PIL not installed, so iteration history will not be generated")
        return
            
        
    # If we are here, then we have imageio and PIL and can make a gif!
        
    # First load the raw images
    raw_images = [file for file in os.listdir(workdir) if (movie_type in file and file.endswith(".png"))]

    # Sort the images based on iteration
    raw_images = sorted(raw_images,key=lambda image: int(image.lstrip(movie_type + "_frame").rstrip('.png')))

    # Create a list to store the imageio objects
    images = []
        
    # Loop through all of the images to draw their iteration number on the plot
    for raw_image in raw_images:
            
        # First get the iteration number
        Iteration_no = raw_image.lstrip("frame").split("_")[0]
            
        # Second, open the raw image
        img = Image.open(PEATSAorder.images_dir + raw_image).convert('RGBA')
            
        # Create the draw object on the image
        draw = ImageDraw.Draw(img,'RGBA')
            
        try:
            # Get the default matplotlib font name
            font_name = findfont(FontProperties(family=[matplotlib.rcParams['font.family'][0]]))
            
            # Setup the matplotlib default font
            fnt = ImageFont.truetype(font_name,20)
        except:
            logging.info("Cant find matplotlib default font, using default and probably small font for gif iteration tag")
                
            # Set up the default font
            fnt = ImageFont.load_default()
            
        # Draw the iteration number
        draw.text((30,10),"Iteration " + Iteration_no,"Black",font=fnt)
            
        # Save the revised image
        img.save(PEATSAorder.images_dir + raw_image.rstrip("png").rstrip(".") + "_with_text.png")
            
        # Add the revised image to the list for imageio to work with
        images.append(imageio.imread(PEATSAorder.images_dir + raw_image.rstrip("png").rstrip(".") + "_with_text.png"))

    # Set the duration of each frame
    kargs = {"duration":1}
        
    # Create the gif
    imageio.mimsave(PEATSAorder.images_dir + movie_type + "_Iteration_history_" + title_string + ".gif",images,"GIF",**kargs)

    # Clean up the working directory
    raw_images = [file for file in os.listdir(PEATSAorder.images_dir) if file.endswith("_with_text.png")]
        
    # Loop through and delete all of the revised images
    for img in raw_images:
        os.remove(PEATSAorder.images_dir + img)

import NLPframe

import os

workdir = 'C:/emtg/EMTG_v9_results/CAESAR_HighFidelity_Sharp100_9propmargin_1em5_RADEC_7262018_185740/'

for filename in os.listdir(workdir):
    if 'NLP_frame' in filename and filename.endswith('.csv'):
        myFrame = NLPframe.NLPframe(workdir + filename)
        myFrame.plot_X()
        myFrame.plot_F()

#make_movie('X')
#make_movie('F')


