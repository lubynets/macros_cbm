alias gsi='ssh -X lubynets@lxi001.gsi.de'
alias mntgsi='sshfs lubynets@lxi001.gsi.de:/u/lubynets/ /home/oleksii/Mount/gsi/'
alias mntlustre2steps='sshfs lubynets@lxi001.gsi.de:/u/lubynets/Mount/Lustre/ /home/oleksii/Mount/Lustre/'
alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@lustre.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
# alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@vae22.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
# alias mntlustre="sshfs -o ssh_command='ssh -J lubynets@lxi001.gsi.de' lubynets@virgo-centos7.hpc.gsi.de:/lustre/cbm/users/lubynets/ /home/oleksii/Mount/Lustre/"
alias mmg="mntgsi && mntlustre && gsi"
alias umntgsi='sudo umount ~/Mount/gsi'
alias umntlustre='sudo umount ~/Mount/Lustre'
alias kssh='killall -9 ssh'

alias cern='ssh -X ollubyne@lxplus.cern.ch'
alias mntcern='sshfs ollubyne@lxplus.cern.ch:/afs/cern.ch/user/o/ollubyne/ /home/user/Mount/cern/'

alias mji='make -j3 install'
alias mj='make -j3'
alias mi='make install'
alias m='make'
alias cd1='cd ..'
alias cd2='cd ../..'
alias cd3='cd ../../..'
alias cd4='cd ../../../..'
alias cd5='cd ../../../../..'
alias cd6='cd ../../../../../..'
alias cd7='cd ../../../../../../..'

alias gst='git status'
alias gpl='git pull'
alias glo='git log'
alias gdf='git diff'
alias gch='git checkout'
alias grv='git remote -v'
alias kt='killall -9 telegram'
alias clang-format='~/.clang-format'
alias pdf2png='pdftoppm -png -cropbox -singlefile'
alias pdf2pngall='pdftoppm -png -cropbox'
alias cpdf='~/.cpdf' # removes pictures from pdf: cpdf -draft input.pdf -o output.pdf

alias AT='source /home/oleksii/soft/AnalysisTree/install_master/bin/AnalysisTreeConfig.sh'
alias qntoolsconfig='source /home/oleksii/soft/.qntoolsconfig.sh'
alias flowdrawconfig='source /home/oleksii/cbmdir/flow_drawing_tools/.config.sh'
alias flowcalcconfig='source /home/oleksii/cbmdir/flow_calculator/.config.sh'
alias qndiscrconfig='source /home/oleksii/cbmdir/qn_discriminator/.config.sh'

alias files644='find ./ -type f -print0 | xargs -0 chmod 644'
alias dirs755='find ./ -type d -print0 | xargs -0 chmod 755'

alias lx='pdflatex main'
alias bx='bibtex main'
alias latexclean='rm main.bbl main.aux main.log main.toc main.pdf main.out main.lot main.lof main.blg'
alias lbll='lx && bx && lx && lx'

alias rootl='root -l'
alias rootlq='root -l -q'
alias rootlb='root -l -b'
alias rootlbq='root -l -b -q'
alias r='root'
alias rl='rootl'
alias rlq='rootlq'
alias rlb='rootlb'
alias rlbq='rootlbq'

########################################################################################################

alias virgo='ssh -Y vae23.hpc.gsi.de'
# alias virgo='virgodebian10'
# alias virgodebian10='ssh -Y virgo-debian10.hpc.gsi.de'
# alias virgocentos7='ssh -Y virgo-centos7.hpc.gsi.de'
alias nano="/lustre/cbm/users/lubynets/soft/nano/install/bin/nano"
# alias mntlustre='sshfs vae22.hpc.gsi.de:/lustre/cbm/users/lubynets/ /u/lubynets/Mount/Lustre'
alias mntlustre='sshfs virgo-centos7.hpc.gsi.de:/lustre/cbm/users/lubynets/ /u/lubynets/Mount/Lustre'
alias cmake='/lustre/cbm/users/lubynets/soft/CMake/install_3.21/bin/cmake'
alias lub='cd /lustre/cbm/users/lubynets'
alias lu='lub'
alias mji='make -j24 install'
alias mj='make -j24'
alias sq='squeue -u lubynets'
alias wsq='watch squeue -u lubynets'
alias sqsc='scancel -u lubynets'
alias cd1='cd ..'
alias cd2='cd ../..'
alias cd3='cd ../../..'
alias cd4='cd ../../../..'
alias cd5='cd ../../../../..'
alias cd6='cd ../../../../../..'
alias cd7='cd ../../../../../../..'
alias cd8='cd ../../../../../../../..'
alias cd9='cd ../../../../../../../../..'
alias gst='git status'
alias gpl='git pull'
alias glo='git log'
alias gdf='git diff'
alias gch='git checkout'
alias tm='tmux'
alias ta='tmux attach'
alias ta0='tmux attach -t 0'
alias ta1='tmux attach -t 1'
alias ta2='tmux attach -t 2'
alias rmlw='rm -rf log/ workdir/'
alias rmwl='rm -rf log/ workdir/'
alias qnanalysisconfig='source /lustre/cbm/users/lubynets/soft/QnAnalysis/install/bin/QnAnalysisConfig.sh'

alias alice='/cvmfs/alice.cern.ch/bin/alienv enter AliPhysics::vAN-20201201_ROOT6-1'
alias checkspace='lfs quota -h -u $USER /lustre/cbm/users/$USER/'
alias countfiles='/lustre/cbm/users/lubynets/.counter_script.sh'

alias rootl='root -l'
alias rootlq='root -l -q'
alias r='root'
alias rl='rootl'
alias rlq='rootlq'
