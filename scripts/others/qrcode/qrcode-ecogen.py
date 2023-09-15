# Create a QR code for the ECOGEN website

# Highly inspired from https://www.geeksforgeeks.org/how-to-generate-qr-codes-with-a-custom-logo-using-python/

# Install qrcode library in virtual environment using the following
# pip3 -m venv venv
# source venv/bin/activate
# pip install -r requirements.txt

# Librairies
import qrcode
from PIL import Image

# Inputs
img = Image.open('logo.png')
url = 'https://code-mphi.github.io/ECOGEN/'

# Taking base width
basewidth = 100
 
# Adjust image size
wpercent = (basewidth/float(img.size[0]))
hsize = int((float(img.size[1])*float(wpercent)))
logo = img.resize((basewidth, hsize), Image.ANTIALIAS)
QRcode = qrcode.QRCode(
    error_correction=qrcode.constants.ERROR_CORRECT_H
)
 
# Adding URL to QRcode
QRcode.add_data(url)
 
# Generating QR code
QRcode.make()
 
# Set colors pixel
QRcolor_pixel = 'Black'
# QRcolor_pixel = (74, 147, 185) # RGB blue
QRcolor_background = 'white'

 
# Adding color to QR code
QRimg = QRcode.make_image(
    fill_color=QRcolor_pixel, back_color=QRcolor_background).convert('RGB')
 
# Set size of QR code + overlay logo (careful img should 
# be less than 30% of QR code to be valid)
pos = ((QRimg.size[0] - logo.size[0]) // 2,
       (QRimg.size[1] - logo.size[1]) // 2)
QRimg.paste(logo, pos)

# Save image
QRimg.save('qrcode.png')
print('QR code generated.')

# Simple version without logo
# qr = qrcode.make(url)
# qr.save('qrcode-ecogen.png')
