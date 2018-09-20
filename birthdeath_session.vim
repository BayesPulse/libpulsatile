let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/Software/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +0 include/bpmod_singlesubject/ss_draw_locations.h
badd +0 src/singlesubject.cpp
badd +0 include/bpmod_singlesubject/birthdeath.h
badd +151 include/bp_datastructures/patient.h
badd +62 include/bp_datastructures/pulseestimates.h
badd +0 include/bp_mcmc/mh.h
badd +678 ~/Projects/BayesPulse/Software/pulsatile/src/mcmc.c
badd +0 src/tempmain.cpp
badd +56 include/bpmod_singlesubject/ss_draw_fixedeffects.h
argglobal
silent! argdel *
$argadd include/bpmod_singlesubject/ss_draw_locations.h
set stal=2
edit include/bpmod_singlesubject/ss_draw_locations.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
3wincmd h
wincmd w
wincmd w
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 31 + 191) / 382)
exe 'vert 2resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 3resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 4resize ' . ((&columns * 116 + 191) / 382)
argglobal
enew
file NERD_tree_1
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal nofen
wincmd w
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 53 - ((52 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
53
normal! 0
wincmd w
argglobal
if bufexists('include/bpmod_singlesubject/birthdeath.h') | buffer include/bpmod_singlesubject/birthdeath.h | else | edit include/bpmod_singlesubject/birthdeath.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 48 - ((47 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
48
normal! 0
wincmd w
argglobal
if bufexists('src/singlesubject.cpp') | buffer src/singlesubject.cpp | else | edit src/singlesubject.cpp | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 159 - ((73 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
159
normal! 05|
wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 191) / 382)
exe 'vert 2resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 3resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 4resize ' . ((&columns * 116 + 191) / 382)
tabedit include/bp_mcmc/mh.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
3wincmd h
wincmd w
wincmd w
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 31 + 191) / 382)
exe 'vert 2resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 3resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 4resize ' . ((&columns * 116 + 191) / 382)
argglobal
enew
file NERD_tree_2
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal nofen
wincmd w
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 72 - ((71 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
72
normal! 029|
wincmd w
argglobal
if bufexists('include/bpmod_singlesubject/ss_draw_locations.h') | buffer include/bpmod_singlesubject/ss_draw_locations.h | else | edit include/bpmod_singlesubject/ss_draw_locations.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 117 - ((73 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
117
normal! 015|
wincmd w
argglobal
if bufexists('~/Projects/BayesPulse/Software/pulsatile/src/mcmc.c') | buffer ~/Projects/BayesPulse/Software/pulsatile/src/mcmc.c | else | edit ~/Projects/BayesPulse/Software/pulsatile/src/mcmc.c | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 513 - ((47 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
513
normal! 032|
wincmd w
3wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 191) / 382)
exe 'vert 2resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 3resize ' . ((&columns * 116 + 191) / 382)
exe 'vert 4resize ' . ((&columns * 116 + 191) / 382)
tabedit src/tempmain.cpp
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 190 + 191) / 382)
exe 'vert 2resize ' . ((&columns * 191 + 191) / 382)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 85 - ((84 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
85
normal! 064|
wincmd w
argglobal
if bufexists('src/singlesubject.cpp') | buffer src/singlesubject.cpp | else | edit src/singlesubject.cpp | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 154 - ((68 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
154
normal! 05|
wincmd w
exe 'vert 1resize ' . ((&columns * 190 + 191) / 382)
exe 'vert 2resize ' . ((&columns * 191 + 191) / 382)
tabedit include/bpmod_singlesubject/birthdeath.h
set splitbelow splitright
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 118 - ((96 * winheight(0) + 51) / 102)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
118
normal! 0
tabnext 2
set stal=1
if exists('s:wipebuf') && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 winminheight=1 winminwidth=1 shortmess=filnxtToO
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
