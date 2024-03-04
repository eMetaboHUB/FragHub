import os

with open(r"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\ALL_DB\MSP\POS\POS_LC.msp", "r") as buffer:
    content = buffer.read()


content = content.split("\n\n")
print(len(content))

content_1 = content[:500000]
content_2 = content[500000:1000000]
content_3 = content[1000000:1500000]
content_4 = content[1500000:2000000]
content_5 = content[2000000:]

content_1 = "\n\n".join(content_1)
content_2 = "\n\n".join(content_2)
content_3 = "\n\n".join(content_3)
content_4 = "\n\n".join(content_4)
content_5 = "\n\n".join(content_5)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part1.msp", "w") as buffer:
    buffer.write(content_1)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part2.msp", "w") as buffer:
    buffer.write(content_2)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part3.msp", "w") as buffer:
    buffer.write(content_3)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part4.msp", "w") as buffer:
    buffer.write(content_4)

with open(r"C:\Users\Axel\Desktop\POS_LC_little_part5.msp", "w") as buffer:
    buffer.write(content_5)

