let SessionLoad = 1
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
cd ~/Projects/BayesPulse/Software/libpulsatile
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +125 tests/mh_tests.cpp
badd +0 include/bp_mcmc/mh.h
badd +140 include/bpmod_singlesubject/birthdeath.h
badd +21 include/bp_mcmc/counter.h
badd +0 include/bp_mcmc/proposalvariance.h
argglobal
silent! argdel *
$argadd tests/mh_tests.cpp
edit tests/mh_tests.cpp
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
exe 'vert 1resize ' . ((&columns * 31 + 159) / 318)
exe 'vert 2resize ' . ((&columns * 95 + 159) / 318)
exe 'vert 3resize ' . ((&columns * 94 + 159) / 318)
exe 'vert 4resize ' . ((&columns * 95 + 159) / 318)
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
let s:l = 119 - ((9 * winheight(0) + 40) / 81)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
119
normal! 0
wincmd w
argglobal
if bufexists('include/bp_mcmc/mh.h') | buffer include/bp_mcmc/mh.h | else | edit include/bp_mcmc/mh.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 85 - ((44 * winheight(0) + 40) / 81)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
85
normal! 05|
wincmd w
argglobal
if bufexists('include/bp_mcmc/proposalvariance.h') | buffer include/bp_mcmc/proposalvariance.h | else | edit include/bp_mcmc/proposalvariance.h | endif
setlocal fdm=syntax
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=2
setlocal fml=1
setlocal fdn=1
setlocal fen
let s:l = 92 - ((11 * winheight(0) + 40) / 81)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
92
normal! 011|
wincmd w
2wincmd w
exe 'vert 1resize ' . ((&columns * 31 + 159) / 318)
exe 'vert 2resize ' . ((&columns * 95 + 159) / 318)
exe 'vert 3resize ' . ((&columns * 94 + 159) / 318)
exe 'vert 4resize ' . ((&columns * 95 + 159) / 318)
tabnext 1
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
