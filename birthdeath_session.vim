let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +150 src/singlesubject.cpp
badd +10 include/counter.h
badd +0 term://.//15141:/bin/bash
badd +0 include/patient.h
badd +54 src/tempmain.cpp
badd +375 include/datastructures.h
badd +0 include/singlesubject/ss_draw_randomeffects.h
badd +8 include/singlesubject/ss_draw_randomeffects_width.h
badd +0 include/singlesubject/ss_draw_fixedeffects.h
argglobal
silent! argdel *
$argadd src/singlesubject.cpp
set stal=2
edit include/singlesubject/ss_draw_randomeffects.h
set splitbelow splitright
wincmd _ | wincmd |
vsplit
wincmd _ | wincmd |
vsplit
2wincmd h
wincmd w
wincmd w
wincmd t
set winminheight=1 winminwidth=1 winheight=1 winwidth=1
exe 'vert 1resize ' . ((&columns * 79 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 79 + 119) / 238)
exe 'vert 3resize ' . ((&columns * 78 + 119) / 238)
argglobal
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 97 - ((29 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
97
normal! 037|
wincmd w
argglobal
if bufexists('include/singlesubject/ss_draw_fixedeffects.h') | buffer include/singlesubject/ss_draw_fixedeffects.h | else | edit include/singlesubject/ss_draw_fixedeffects.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 74 - ((24 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
74
normal! 07|
wincmd w
argglobal
if bufexists('term://.//15141:/bin/bash') | buffer term://.//15141:/bin/bash | else | edit term://.//15141:/bin/bash | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 988 - ((50 * winheight(0) + 28) / 56)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
988
normal! 03|
wincmd w
exe 'vert 1resize ' . ((&columns * 79 + 119) / 238)
exe 'vert 2resize ' . ((&columns * 79 + 119) / 238)
exe 'vert 3resize ' . ((&columns * 78 + 119) / 238)
tabedit include/patient.h
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
let s:l = 1 - ((0 * winheight(0) + 46) / 93)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
tabnext 1
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
