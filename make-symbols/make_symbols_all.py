import os
import string
from PIL import Image, ImageDraw, ImageFont

# Path to your installed (or copied) Agathodaimon font:
FONT_PATH = "/usr/share/fonts/truetype/agathodaimon/agathodaimon.ttf"

def generate_alphabet_images(output_dir="alphabet_images"):
    """
    Generates 500x500 PNGs of A-Z, a-z using the Agathodaimon font.
    """
    # Make sure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # We'll take uppercase A-Z + lowercase a-z
    # from Python's built-in 'string' module
    alphabet = string.ascii_uppercase + string.ascii_lowercase

    # For each character, create an image and save it
    for ch in alphabet:
        # Create a 500x500 RGBA image with a transparent background
        img = Image.new("RGBA", (500, 500), (255, 255, 255, 0))
        draw = ImageDraw.Draw(img)

        # Load the font at a size that fits well within 500x500
        font = ImageFont.truetype(FONT_PATH, 400)

        # Measure the text size (width, height) so we can center it
        w, h = draw.textsize(ch, font=font)
        x = (500 - w) // 2
        y = (500 - h) // 2

        # Draw the letter in black (0,0,0,255)
        draw.text((x, y), ch, font=font, fill=(0, 0, 0, 255))

        # Save the image to disk (letter as filename)
        out_filename = os.path.join(output_dir, f"{ch}.png")
        img.save(out_filename)

    print(f"Images saved to: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    generate_alphabet_images()

