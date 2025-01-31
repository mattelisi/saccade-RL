from PIL import Image, ImageDraw, ImageFont

# Path to your Agathodaimon font file
FONT_PATH = "/usr/share/fonts/truetype/agathodaimon/agathodaimon.ttf"

def render_symbol(symbol, out_filename, font_size=400, img_size=500):
    # Create a blank RGBA image (transparent background)
    img = Image.new("RGBA", (img_size, img_size), (255, 255, 255, 0))
    draw = ImageDraw.Draw(img)

    # Load the font
    font = ImageFont.truetype(FONT_PATH, font_size)

    # Measure text size to center it
    text_width, text_height = draw.textsize(symbol, font=font)

    # Coordinates so that text is centered
    x = (img_size - text_width) / 2
    y = (img_size - text_height) / 2

    # Draw the symbol in black (0,0,0,255)
    draw.text((x, y), symbol, font=font, fill=(0, 0, 0, 255))

    # Save as PNG (with alpha channel intact)
    img.save(out_filename)

# Usage example
if __name__ == "__main__":
    render_symbol("A", "agatho_A.png")
    render_symbol("B", "agatho_B.png")
