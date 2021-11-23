from PIL import Image

img = Image.open('Y:/common/_Personal_Folder/XS/esi_spray.png')
img = img.convert("RGBA")
datas = img.getdata()

newData = []
for item in datas:
    # if item[0] == 0 and item[1] == 0 and item[2] == 0:
    if item[0] == item[1] == item[2]:
        newData.append((255, 255, 255, 0))
    else:
        newData.append(item)

img.putdata(newData)
img.save("Y:/common/_Personal_Folder/XS/esi_spray_transparent.png", "PNG")

