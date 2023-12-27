# 1. 滞留在Cutadapt步骤
停止在cutadapt步骤，等候一天无任何反应。
逐行运行代码发现，原因是停留在fastqc一步，尤其是对于较大的样本（目前发现的是双端150bp，文件大小各自为2.9G），报错信息：

>>> Now running FastQC on the validated data BratOE-Rep2_val_1.fq.gz<<<

Started analysis of BratOE-Rep2_val_1.fq.gz
Exception in thread "Thread-1" java.lang.OutOfMemoryError: Java heap space
	at uk.ac.babraham.FastQC.Utilities.QualityCount.<init>(QualityCount.java:33)
	at uk.ac.babraham.FastQC.Modules.PerTileQualityScores.processSequence(PerTileQualityScores.java:281)
	at uk.ac.babraham.FastQC.Analysis.AnalysisRunner.run(AnalysisRunner.java:89)
	at java.lang.Thread.run(Thread.java:745)

经查发现：https://github.com/s-andrews/FastQC/issues/86
使用代码 'fastqc BratOE-Rep2_val_1.fq.gz -t 24'错误会消失！
总结： 使用流程时，注意双端测序数据的每个文件大小不要大于2.8G！

