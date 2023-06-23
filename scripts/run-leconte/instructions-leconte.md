```console
kol@leconte:~$ git clone https://github.com/JuliaORNL/GrayScott.jl.git
kol@leconte:~$ cp -r GrayScott.jl/scripts/run-leconte/run .
kol@leconte:~$ cd run
kol@leconte:~/run$ chmod +x config_leconte.sh job_leconte.sh
kol@leconte:~/run$ source config_leconte.sh
kol@leconte:~/run$ # read job script documentation
kol@leconte:~/run$ cat job_leconte.sh
kol@leconte:~/run$ mkdir 001
kol@leconte:~/run$ cp $GS_DIR/examples/settings-files.json 001
kol@leconte:~/run$ ./job_leconte.sh 001
```
