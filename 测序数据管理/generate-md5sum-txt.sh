#!/bin/bash
# Using getopt

# 功能：
# 1. 计算目录下（-d）所有文件的md5sum值（仅仅计算不存在MD5.txt以及md5sum值的样本），
#    并保存到同一目录下，文件名为md5sum.txt(-n)(每行：文件名\tabMd5sum值)
#
# 2. check 功能（-c）(if set true, -f will not work!)
# 2.1. 列出哪些文件的md5sum值不存在。
# 2.2. 列出MD5.txt中不存在的文件。
#
# 3. -f:true 强制重新计算md5sum值，并保存，如1。
# 4. 文件后缀 -t:fastq.gz 

# 注意
# 1. 目录不允许为软连接，请进入到非链接后的目录，比如cd  /data5/Xilab-Lab-Data-Backup/ChenJun会正常运行，
# 但是cd /home/xilab/Data_Backup/ChenJun中再去运行命令是错误的！因为这里的ChenJun是一个软连接！

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e
# Default parameters
target_dir=`pwd`
output_name='MD5.txt'
input_pattern='fastq.gz'
check_md5sum=false
force_recalculate=false

# Augument Parsing
print_usage_and_exit(){
    echo "Usage: $0 
            [-d <working directory> = $target_dir default]
            [-n <output file name> = 'MD5.txt' default]
            [-t <input sample postfix> = 'fastq.gz' default]
            [-c <list files without md5sum value> = false default] 
            [-f <force recalculate md5sum value> = false default, No effects when run with '-c = true']
            [-h <print help info>]
        Use cases:
            # check samples and md5 files: 
            sh generate-md5sum-txt.sh -c true
            # first time use or just want to regenerate all
            sh generate-md5sum-txt.sh -f true
            # just run on newlly added files: 
            sh generate-md5sum-txt.sh
            "
    exit 1
}

while getopts ":d:n:t:c:f:" opt; do
    case $opt in
        d)
            target_dir="$OPTARG"
            echo "-d <working directory> = $target_dir"
            ;;
        n)
            output_name="$OPTARG"
            echo "-n <output file name> = $output_name"
            ;;
        t)
            input_pattern="$OPTARG"
            echo "-t <input sample postfix> = $input_pattern"
            ;;
        c)
            check_md5sum="$OPTARG"
            echo "-c <list files without md5sum value> = $check_md5sum"
            ;;
        f)
            force_recalculate="$OPTARG"
            echo "<force recalculate md5sum value> = $force_recalculate"
            ;;
        h)
            print_usage_and_exit
            ;;
        \?)
            echo "Invalid option: -$OPTARG"
            print_usage_and_exit
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            print_usage_and_exit
            ;;
    esac
done

echo ">>> Working Directory: ${target_dir}..."
echo ">>> File pattern is set to ${input_pattern}..."

# 列出目录下的所有符合pattern的文件，并保存到files.path文件中。
## 检查是否有已经存在的files.path文件
if [[ -f $target_dir/files.path ]]; then
    echo "Found 'files.path' file under the direcotry: $target_dir"
    while true; do
        read -p "Do you wish to remove this file?" yn
        case $yn in
            [Yy]* ) rm $target_dir/files.path; break;;
            [Nn]* ) echo "Please check the file - files.path."; exit 1;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi
## 检查是否有已经存在的MD5.path文件
if [[ -f $target_dir/MD5.path ]]; then
    echo "Found 'MD5.path' file under the direcotry: $target_dir"
    while true; do
        read -p "Do you wish to remove this file?" yn
        case $yn in
            [Yy]* ) rm $target_dir/MD5.path; break;;
            [Nn]* ) echo "Please check the file - MD5.path."; exit 1;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi
## 查找所有符合pattern的文件，并保存到files.path和MD5.path
### files.path
find $target_dir -name *${input_pattern} >> ${target_dir}/files.path
samples_numbers=`cat $target_dir/files.path | wc -l`
if [[ $samples_numbers -ge 1 ]]; then
    echo ">>> $samples_numbers samples found..."
else
    echo ">>> zero samples found by using pattern - ${input_pattern}, please check the file pattern used..."
    exit 1
fi
### MD5.path
find $target_dir -name $output_name >> ${target_dir}/MD5.path
md5txt_numbers=`cat $target_dir/MD5.path | wc -l`
if [[ $md5txt_numbers -ge 1 ]]; then
    echo ">>> $md5txt_numbers md5 files found..."
else
    echo ">>> No $output_name files found, you may need to check the file pattern used..."
fi


# functions
if [[ "${check_md5sum}" = true ]]; then
    echo ">>> Start to check..."
    ## check function
    ### file level
    echo ">>> checking at file level, files without md5 value will be show bellow..."
    cat ${target_dir}/files.path | while read line
    do
        fq_dir=`dirname $line`
        fq_name=`basename $line`
        # 1. 检查同目录下MD5.txt是否存在，存在的话，检查MD5.txt文件的第三列是否有该文件名，如果没有则输出该文件（含路径）
        if [[ -f ${fq_dir}/${output_name} ]]; then
            echo $(cut -f3 -d' ' ${fq_dir}/${output_name}) | grep -wq "$fq_name" &>/dev/null || echo $line
            # grep -xq为绝对匹配，不适用于有多个元素的array
            # 使用-wq的可能问题，基于正则匹配，非绝对匹配，可能会有潜在错误，比如后期修改了fastq的文件名。
        else
            echo $line
        fi
    done
    ### md5.txt level
    echo ">>> checking at md5.txt level, md5.txt without existed files will be show bellow..."
    cat ${target_dir}/MD5.path | while read line
    do
        md5_dir=`dirname $line`
        for i in $(cut -f3 -d' ' $line)
        do
            [[ -f ${md5_dir}/${i} ]] || echo "$line - ${i}"
        done
    done
    echo "If you see any files listed above, you have to check the md5 file manually to delete the extra enntry."
else
    if [[ "$force_recalculate" = true ]]; then
        echo ">>> Running as force recalculate mode..."
        pre_wk_dir=""
        cat ${target_dir}/files.path | while read line
        do
            echo "Processing sample $line"
            wk_dir=`dirname $line`
            if [[ "$wk_dir" != "$pre_wk_dir" ]]; then
                cd $wk_dir
                if [[ -f ${output_name} ]]; then
                    rm  ${output_name}
                fi
            fi
            md5sum `basename $line` >> ${output_name}
            pre_wk_dir=${wk_dir}
        done
    else
        ## only handle the newly added files
        echo ">>> Running as newly added mode..."
        cat ${target_dir}/files.path | while read line
        do
            fq_dir=`dirname $line`
            fq_name=`basename $line`
            # 1. 检查同目录下MD5.txt是否存在，存在的话，检查MD5.txt文件的第三列是否有该文件名，如果没有则输出该文件（含路径）
            if [[ -f ${fq_dir}/${output_name} ]]; then
                echo $(cut -f3 -d' ' ${fq_dir}/${output_name}) | grep -wq "$fq_name" &>/dev/null || echo $line >> ${target_dir}/tmp.run.generate.md5sum.list
            else
                echo $line >> ${target_dir}/tmp.run.generate.md5sum.list
            fi
        done
        if [[ ! -f ${target_dir}/tmp.run.generate.md5sum.list ]]; then
            echo ">>> Now newly added files found, Exit..."
            exit 1
        fi
        cat ${target_dir}/tmp.run.generate.md5sum.list | while read new
        do
            echo "Processing sample: ${new}"
            cd `dirname ${new}`
            md5sum `basename ${new}` >> ${output_name}
        done
        rm ${target_dir}/tmp.run.generate.md5sum.list
        # <<COMMENT
        #     Not tested!
        #     new_fq=()
        #     for Afile in `cat ${target_dir}/files.path`
        #     do
        #         fq_dir=`dirname $Afile`
        #         fq_name=`basename $Afile`
        #         # 1. 检查同目录下MD5.txt是否存在，存在的话，检查MD5.txt文件的第三列是否有该文件名，如果没有则输出该文件（含路径）
        #         if [[ -f ${fq_dir}/${output_name} ]]; then
        #             echo $(cut -f3 -d' ' ${fq_dir}/${output_name}) | grep -wq "$fq_name" &>/dev/null || new_fq=("${new_fq[@]}" $Afile)
        #         else
        #             new_fq=("${new_fq[@]}" $Afile)
        #         fi
        #     done
        #     echo ${new_fq[@]}
        #     for Anew_fq in $new_fq
        #     do
        #         echo "Processing sample $Anew_fq"
        #         cd `dirname $Anew_fq`
        #         md5sum `basename $Anew_fq` >> ${output_name}
        #     done
        # COMMENT
    fi    
fi   
# <<COMMENT
#     有个问题，就是处理新加样本时，定义一个new_fq=()
#     然后在while循环中，把新加入的文件添加到这个变量中 # new_fq=("${new_fq[@]}" $line)
#     问题是，推出循环后，new_fq()的值还是空的！ 咋回事呢?
#     之所以我写的脚本出现了输出是空的问题，原因就在这里
#     linux执行shell时，会创建“子shell”运行shell中的命令，当运行到非内建指令时，会创建“孙shell”运行非内建指令
#     变量的作用于在每个shell中有效，所以，非内建指令中定义的这些变量就只能在孙shell运行，而在子shell中不生效，所以，即便我在while中给path_all赋值了，子shell中也不会获取到这个值。
#     解决这个问题的办法有两种，如下
#     如果不是必须使用管道符的方式写while循环，可以用重定向的写法，这种写法循环内的变量在子shell中是生效的，比较简便
#     如果非要使用管道符的方式，可以创建临时文件，用于存放孙shell中的输出
#     试一下for循环？
# COMMENT

rm ${target_dir}/MD5.path
rm ${target_dir}/files.path
trap : 0