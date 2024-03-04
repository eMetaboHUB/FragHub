

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3_4_3_4_5_1_3_4.msp", "r") as buffer:
    content = buffer.read()

content = content.split("\n\n")
print(len(content))

content_1 = content[0]
content_2 = content[1]
content_3 = content[2]
content_4 = content[3]
content_5 = content[4]

# content_1 = "\n\n".join(content_1)
# content_2 = "\n\n".join(content_2)
# content_3 = "\n\n".join(content_3)
# content_4 = "\n\n".join(content_4)
# content_5 = "\n\n".join(content_5)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3_4_3_4_5_1_3_4_1.msp", "w") as buffer:
    buffer.write(content_1)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3_4_3_4_5_1_3_4_2.msp", "w") as buffer:
    buffer.write(content_2)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3_4_3_4_5_1_3_4_3.msp", "w") as buffer:
    buffer.write(content_3)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3_4_3_4_5_1_3_4_4.msp", "w") as buffer:
    buffer.write(content_4)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3_4_3_4_5_1_3_4_5.msp", "w") as buffer:
    buffer.write(content_5)
