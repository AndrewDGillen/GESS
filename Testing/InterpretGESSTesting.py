import glob
import csv
from statistics import stdev, mean
from collections import defaultdict

def main():

    testresults = []

    for file in glob.glob('/home/ag_jdowlab/Desktop/FlyGene/GESS_BATCH_OUTPUTS/*.log'):
        
        try:
            source, gesstype, dtype, _ = file.split('/')[-1].split('_')
        except:
            continue
        dtypecount = int(dtype)
        
        with open(file) as testdata_in:

            linecount = 0
            runningtotal = defaultdict(list)

            for line in testdata_in:
                linecount += 1

                if linecount % 10 in [7, 8, 9]:
                    timetype, time = line.strip().split('\t')
                    
                    time_split = time[:-1].split('m')
                    actual_time = (float(time_split[0])*60) + (float(time_split[1]))

                    runningtotal[timetype].append(actual_time)

            indiv_test_result = [source, gesstype, dtypecount]

            for timetype in ['real', 'user', 'sys']:
                avg_runtime = mean(runningtotal[timetype])
                stdev_runtime = stdev(runningtotal[timetype])

                indiv_test_result.extend([avg_runtime, stdev_runtime])

        
        testresults.append(indiv_test_result)

    testresults=sorted(testresults,key=lambda x:(x[0], x[1], x[2]))
    
    with open('/home/ag_jdowlab/Desktop/FlyGene/GESS_UNITTESTING.csv', 'w+') as unittestfile:
        OUT = csv.writer(unittestfile)

        OUT.writerow(['Data source','GESS type', 'Categories', 'REAL Runtime', 'REAL StDev', 'USER runtime', 'USER StDev', 'SYS runtime', 'SYS StDev'])

        for outline in testresults:
            OUT.writerow(outline)
            
if __name__ == "__main__":
    main()