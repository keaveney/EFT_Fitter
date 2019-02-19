######insert "preliminary" tag in .eps file
import os
import subprocess

path_1 = "/Users/keaveney/docs/TOP-17-014-PAS/myDir/notes/TOP-17-014/trunk/fig_logo/diff_results/preliminary/normalized/particle/"
path_2 = "/Users/keaveney/docs/TOP-17-014-PAS/myDir/notes/TOP-17-014/trunk/fig_logo/diff_results/preliminary/normalized/parton/"
path_3 = "/Users/keaveney/docs/TOP-17-014-PAS/myDir/notes/TOP-17-014/trunk/fig_logo/diff_results/preliminary/absolute/particle/"
path_4 = "/Users/keaveney/docs/TOP-17-014-PAS/myDir/notes/TOP-17-014/trunk/fig_logo/diff_results/preliminary/absolute/parton/"

path_list = [path_1, path_2, path_3, path_4]
searchstring = "gsave  2268 2176 0 0 C 691.54 1986 t 0 r /Helvetica-Bold findfont 91.3174 sf 0 0 m (CMS) show NC gr"
insertstring = "gsave  2268 2176 0 0 C 909.37 1993.76 t 0 r /Helvetica-Oblique findfont 79.9028 sf 0 0 m (Preliminary) show NC gr"


for path in path_list:
    for filename in os.listdir(path):
        if filename.endswith(".eps"):
            targetfile = os.path.join(path, filename)
            print(targetfile)
            inputfile = open(targetfile, 'r').readlines()
            write_file = open(targetfile,'w')
            for line in inputfile:
                write_file.write(line)
                if searchstring in line:
                    write_file.write(insertstring + "\n")
            write_file.close()
            subprocess.call(['epstopdf', targetfile])
            continue
        else:
            continue

