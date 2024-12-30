import os
def file_name(file_dir):
    for root, dirs, files in os.walk(file_dir):
        print(root)
        print(dirs)
        print(files)

file_name('C:\Users\Xi_Lab\PycharmProjects\Test')