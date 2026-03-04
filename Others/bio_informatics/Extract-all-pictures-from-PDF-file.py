import fitz   # install - PyMuPDF
import time
import re
import os


def pdf2pic(path, pic_path):
    t0 = time.process_time()  # 生成图片初始时间
    checkXO = r"/Type(?= */XObject)"  # 使用正则表达式来查找图片
    checkIM = r"/Subtype(?= */Image)"
    doc = fitz.open(path)  # 打开pdf文件
    imgcount = 0  # 图片计数
    lenXREF = doc.xref_length()
    # lenXREF = doc._getXrefLength()  # 获取对象数量长度
    #lenXREF = doc._getObjectString()

    # 打印PDF的信息
    print("文件名:{}, 页数: {}, 对象: {}".format(path, len(doc), lenXREF - 1))

    # 遍历每一个对象
    for i in range(1, lenXREF):
        text = doc.xrefObject(i)    # 定义对象字符串
        isXObject = re.search(checkXO, text)  # 使用正则表达式查看是否是对象
        isImage = re.search(checkIM, text)  # 使用正则表达式查看是否是图片
        print(text)
        # if not isXObject or not isImage:  # 如果不是对象也不是图片，则continue
        if not isImage:  # 如果不是对象也不是图片，则continue
            continue
        # print(text)
        imgcount += 1
        pix = fitz.Pixmap(doc, i)  # 生成图像对象
        new_name = "图片{}.png".format(imgcount)  # 生成图片的名称
        if pix.n < 5:  # 如果pix.n<5,可以直接存为PNG
            print("pix.n < 5")
            try:
                pix.writePNG(os.path.join(pic_path, new_name))
            except RuntimeError:  # New code added!
                print("RuntimeError Happened! trying new code!")
                pix0 = fitz.Pixmap(fitz.csRGB, pix)
                pix0.writePNG(os.path.join(pic_path, new_name))
            except ValueError:
                pix0 = fitz.Pixmap(fitz.csRGB, pix)
                pix0.writePNG(os.path.join(pic_path, new_name))
        else:  # 否则先转换CMYK
            print("pix.n >= 5")
            pix0 = fitz.Pixmap(fitz.csRGB, pix)
            pix0.writePNG(os.path.join(pic_path, new_name))
        pix = None  # 释放资源
        t1 = time.process_time()  # 图片完成时间static
        print("运行时间:{}s".format(t1 - t0))
        print("提取了{}张图片".format(imgcount))


def main():
    road = "C:/Users/Nobody/Desktop/干细胞专项/"
    pdf_files = []
    for i in os.listdir(road):
        if os.path.isfile(road + i):
            pdf_files.append(i)
    j = 1
    for i in pdf_files:
        path = road + i
        pic_path = road + 'Extracted-Pictures-From-File-' + os.path.splitext(i)[0]
        if os.path.exists(pic_path):
            os.removedirs(pic_path)
        os.mkdir(pic_path)
        pdf2pic(path, pic_path)
        j = j + 1
        print("\n")


main()