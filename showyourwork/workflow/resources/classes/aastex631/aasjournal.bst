%% AASJournal.bst June 3, 2019
%%
%% Revision 1.18 2019/6/3 Amy Hendrickson
%% Took out `v' preceding version number in @misc
%%
%% Rewrote format.pid to implement this precedence:
%% if doi, print that, if eprint print that, if neither
%% and there is a url, use the url

%% Took out format.url output generally because it is included in format.pid,
%% so undesirable to call for it a second time.

%%
%% Revision 1.17 2019/3/15 Amy Hendrickson
%% Added DOI to every field where format.url was found. Kept format.url in case
%% it is needed.
%%

%% Revision 1.16 2017/10/26 Amy Hendrickson
%% `url' only active when used in misc
%%
%% `note' is disabled in all cases; as is adsnote and adsurl.
%%
%% Repair to @misc call to aaspublisher to prevent `, ,'
%% when there is no publisher listed.
%%
%% Make correct hyperlinked content for doi, ascl, and arXiv;
%% changed definition of FUNCTION {format.eprint} to compare
%%  archivePrefix with "arXiv" and when they match do
%%        ' {new.block "\doarXiv{" eprint * "}" * }
%%                      otherwise do
%%        ' {new.block "\doeprint{" eprint * "}" * }
%%
%% \dodoi, \doarXiv and \doeprint definitions are sent to the top of the .bbl file
%% with commands found in FUNCTION {begin.bib}.
%%
%% bbl.version no longer used, could be deleted.
%%
%% Version number is added when author enters: version = {number}, only in `misc'
%%
%% Removed `format.doi output' in {article} because in 
%% combination with { format.pid output } we could get double doi entries 
%% in the bibliography. 
%%
%% Rewrote format.pid so that we can check to see if we have a doi number; in
%% which case we use that and suppress the arXiv and eprint combination. However,
%% we do Not eliminate the url field if the author has entered it.
%%
%% Deleted format.miscversion and made a simplified way to produce version in @misc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Revision 1.15 2016/10/5 Amy Hendrickson
%% Changed bbl.version so that there is no longer a trailing `.'

%% Revision 1.14 2016/9/14 Amy Hendrickson
%% Added URL field 

%% Revision 1.13 2016/8/9 Amy Hendrickson
%% Change to presort and bibsort to get multiple entries with same author/year
%% combination to be ordered in bibliography in same order as \cite{}s entered in text.

%% Revision 1.12 2016/2/11 Amy Hendrickson
%% Changes to @misc function to specifications by Dr. August Muench.
%% Made comma appear between title and publisher whether or not there is a version number.

%% Revision 1.11 2015/12/31 Amy Hendrickson
%% Changes to @misc function to specifications by Dr. August Muench.
%% Added version and publisher fields.

%% Revision 1.10 2013/03/11 13:01:55 alberto
%% Reworked code to use format.pid for articles without proper
%% pagination which then leads the system to use either an eprint
%% id or a DOI; this will solve most problems with citing early
%% access articles (no volume of page available, but online).
%% Removed bibliography count from preamble since it seems to
%% confuse some versions of latex.
%% 
%% Revision 1.9 2013/2/19 Erik Tollerud <erik.tollerud@yale.edu> 
%% fixes for spurious error messages due to rev 1.8 changes
%% also fixed missing publisher field printing in "book"
%%
%% Revision 1.8  2013/01/10 16:31:34  alberto
%% Added proper formatting of arXiv eprint entries output by ADS
%% (these are article records, but the eprint id is printed only 
%% if no traditional page number is available).  Also added support
%% for ASCL entries under the "misc" category.
%%
%% Revision 1.7  2011/11/21 13:10:12  alberto
%% Implemented new ApJ formatting rules: for papers with more than
%% 5 authors, print first 3 followed by "et al."
%%
%% Revision 1.6  2009/05/18 15:13:02  alberto
%% Fixed sorting of references with truncated lists of authors
%% in sort.format.names
%% Fix courtesy of Stephen Fegan <sfegan@llr.in2p3.fr>
%%
%% Revision 1.5  2009/03/19 13:15:11  alberto
%% Truncate \bibitem entries with more than 8 authors.
%% Fix courtesy of Craig Wiegert <wigie@alum.mit.edu>
%% and Andrew R. Marble <amarble@as.arizona.edu>
%%
%% Revision 1.4  2002/06/18 16:37:48  alberto
%% Add comma after first author in two-author reference
%% Fix courtesy of Tim Robishaw <robishaw@astron.berkeley.edu>
%%
%% Revision 1.3  2000/04/20 22:17:50  jbaker
%% Fixed INBOOK bug, now works essentially like BOOK.
%%
%% Revision 1.2  1998/08/30 22:35:45  jbaker
%% Added RCS keywords.
%%
%% 1998/08/12  J Baker
%% Tweaked by hand to get correct results for ApJ.  Added functions from
%% astrobib.


%%
%% This originally was file `apj.bst',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% merlin.mbs  (with options: `,ay,nat,nm-rev,nmdash,dt-beg,yr-per,note-yr,atit-u,jtit-x,jttl-rm,thtit-a,vnum-x,volp-com,jpg-1,pp-last,btit-rm,add-pub,pub-par,pre-edn,edby,edbyx,blk-com,fin-bare,ppx,ed,abr,ord,jabr,amper,em-x')
%% ----------------------------------------
%% *** Bibliographic Style for ApJ ***
%% 
 %-------------------------------------------------------------------
 % The original source file contains the following version information:
 % \ProvidesFile{merlin.mbs}[1998/02/25 3.85a (PWD)]
 %
 % NOTICE:
 % This file may be used for non-profit purposes.
 % It may not be distributed in exchange for money,
 %   other than distribution costs.
 %
 % The author provides it `as is' and does not guarantee it in any way.
 %
 % Copyright (C) 1994-98 Patrick W. Daly
 %-------------------------------------------------------------------
 %   For use with BibTeX version 0.99a or later
 %-------------------------------------------------------------------
 % This bibliography style file is intended for texts in ENGLISH
 % This is an author-year citation style bibliography. As such, it is
 % non-standard LaTeX, and requires a special package file to function properly.
 % Such a package is    natbib.sty   by Patrick W. Daly
 % The form of the \bibitem entries is
 %   \bibitem[Jones et al.(1990)]{key}...
 %   \bibitem[Jones et al.(1990)Jones, Baker, and Smith]{key}...
 % The essential feature is that the label (the part in brackets) consists
 % of the author names, as they should appear in the citation, with the year
 % in parentheses following. There must be no space before the opening
 % parenthesis!
 % With natbib v5.3, a full list of authors may also follow the year.
 % In natbib.sty, it is possible to define the type of enclosures that is
 % really wanted (brackets or parentheses), but in either case, there must
 % be parentheses in the label.
 % The \cite command functions as follows:
 %   \citet{key} ==>>                Jones et al. (1990)
 %   \citet*{key} ==>>               Jones, Baker, and Smith (1990)
 %   \citep{key} ==>>                (Jones et al., 1990)
 %   \citep*{key} ==>>               (Jones, Baker, and Smith, 1990)
 %   \citep[chap. 2]{key} ==>>       (Jones et al., 1990, chap. 2)
 %   \citep[e.g.][]{key} ==>>        (e.g. Jones et al., 1990)
 %   \citep[e.g.][p. 32]{key} ==>>   (e.g. Jones et al., p. 32)
 %   \citeauthor{key} ==>>           Jones et al.
 %   \citeauthor*{key} ==>>          Jones, Baker, and Smith
 %   \citeyear{key} ==>>             1990
 %---------------------------------------------------------------------

ENTRY
  { address
    archivePrefix
    author
    booktitle
    chapter
    doi
    edition
    editor
    eprint
    howpublished
    institution
    journal
    key
    month
    number
    organization
    pages
    publisher
    school
    series
    title
    misctitle
    type
    volume
    year
    version
    url
  }
  {}
  { label extra.label sort.label short.list }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { add.period$ write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { add.period$ " " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {fin.entry}
{ duplicate$ empty$
    'pop$
    'write$
  if$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
        'skip$
        { after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {add.blank}
{  " " * before.all 'output.state :=
}

FUNCTION {date.block}
{
  skip$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ skip$ }

FUNCTION {capitalize}
{ "u" change.case$ "t" change.case$ }

FUNCTION {space.word}
{ " " swap$ * " " * }

 % Here are the language-specific definitions for explicit words.
 % Each function has a name bbl.xxx where xxx is the English word.
 % The language selected here is ENGLISH
FUNCTION {bbl.and}
{ "and"}

FUNCTION {bbl.editors}
{ "eds." }

FUNCTION {bbl.editor}
{ "ed." }

FUNCTION {bbl.edby}
{ "edited by" }

FUNCTION {bbl.edition}
{ "edn." }

FUNCTION {bbl.volume}
{ "Vol." }

FUNCTION {bbl.version}
{ "" }

FUNCTION {bbl.of}
{ "of" }

FUNCTION {bbl.number}
{ "no." }

FUNCTION {bbl.nr}
{ "no." }

FUNCTION {bbl.in}
{ "in" }

FUNCTION {bbl.pages}
{ "" }

FUNCTION {bbl.page}
{ "" }

FUNCTION {bbl.chapter}
{ "Ch." }

FUNCTION {bbl.techrep}
{ "Tech. Rep." }

FUNCTION {bbl.mthesis}
{ "Master's thesis" }

FUNCTION {bbl.phdthesis}
{ "PhD thesis" }

FUNCTION {bbl.first}
{ "1st" }

FUNCTION {bbl.second}
{ "2nd" }

FUNCTION {bbl.third}
{ "3rd" }

FUNCTION {bbl.fourth}
{ "4th" }

FUNCTION {bbl.fifth}
{ "5th" }

FUNCTION {bbl.st}
{ "st" }

FUNCTION {bbl.nd}
{ "nd" }

FUNCTION {bbl.rd}
{ "rd" }

FUNCTION {bbl.th}
{ "th" }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"Jun."}

MACRO {jul} {"Jul."}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sep."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

FUNCTION {eng.ord}
{ duplicate$ "1" swap$ *
  #-2 #1 substring$ "1" =
     { bbl.th * }
     { duplicate$ #-1 #1 substring$
       duplicate$ "1" =
         { pop$ bbl.st * }
         { duplicate$ "2" =
             { pop$ bbl.nd * }
             { "3" =
                 { bbl.rd * }
                 { bbl.th * }
               if$
             }
           if$
          }
       if$
     }
   if$
}

INTEGERS { nameptr namesleft numnames numtrunc }

FUNCTION {format.names}
{ 's :=
  #0 'numtrunc :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
  numnames #5 > 
    { #3 'numtrunc := }
    {}
  if$
  { namesleft #0 > }
  { s nameptr
    "{vv~}{ll}{, jj}{, f.}" format.name$
    't :=
    nameptr #1 >
      {
	nameptr #1 - numtrunc = 
	  { 
	    #1 'namesleft :=
	    "others" 't :=
	  }
	  {}
	if$

        namesleft #1 >
          { ", " * t * }
          {
            numnames #1 >
              { "," * }
              'skip$
            if$
            s nameptr "{ll}" format.name$ duplicate$ "others" =
              { 't := }
              { pop$ }
            if$
            t "others" =
              {
                " {et~al.}" *
              }
              { " \& " * t * }
            if$
          }
        if$
      }
      't
    if$
    nameptr #1 + 'nameptr :=
    namesleft #1 - 'namesleft :=
  }
  while$
}

FUNCTION {format.names.ed}
{ 's :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{f.~}{vv~}{ll}{, jj}"
      format.name$
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              numnames #2 >
                { "," * }
                'skip$
              if$
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {
                  " {et~al.}" *
                }
                { " \& " * t * }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {format.key}
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
        { ", " * bbl.editors * }
        { ", " * bbl.editor * }
      if$
    }
  if$
}

FUNCTION {format.in.editors}
{ editor empty$
    { "" }
    { editor format.names.ed
    }
  if$
}


FUNCTION {format.title}
{ title empty$
    { "" }
    { title
    }
  if$
}

FUNCTION {format.full.names}
{'s :=
  #1 'nameptr :=
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr
      "{vv~}{ll}" format.name$
      't :=
      nameptr #1 >
        {
          namesleft #1 >
            { ", " * t * }
            {
              numnames #2 >
                { "," * }
                'skip$
              if$
              s nameptr "{ll}" format.name$ duplicate$ "others" =
                { 't := }
                { pop$ }
              if$
              t "others" =
                {
                  " {et~al.}" *
                }
                { " \& " * t * }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {author.editor.key.full}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            'key
          if$
        }
        { editor format.full.names }
      if$
    }
    { author format.full.names }
  if$
}

FUNCTION {author.key.full}
{ author empty$
    { key empty$
         { cite$ #1 #3 substring$ }
          'key
      if$
    }
    { author format.full.names }
  if$
}

FUNCTION {editor.key.full}
{ editor empty$
    { key empty$
         { cite$ #1 #3 substring$ }
          'key
      if$
    }
    { editor format.full.names }
  if$
}

FUNCTION {make.full.names}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.full
    { type$ "proceedings" =
        'editor.key.full
        'author.key.full
      if$
    }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem[{" write$
  label write$
  ")" make.full.names duplicate$ short.list =
     { pop$ }
     { * }
   if$
  "}]{" * write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}


FUNCTION {format.url}
{ url empty$
    { "" }
    { new.block "\url{" url * "}" * }
  if$
}

FUNCTION {format.eprint}
{ eprint empty$
    { "" }
    { archivePrefix empty$
        { "" }
        { archivePrefix "arXiv" =
        {new.block "\doarXiv{" eprint * "}" * }
        {new.block "\doeprint{" eprint * "}" * }
        if$}
      if$
    }
  if$
}

% Formats a persistent identifier (PID)
% (typically a DOI or an eprint id)
% You can easily swap the order of which PID
% gets picked
%

FUNCTION {format.doi}
{ doi empty$
    { "" }
    { "\dodoi{" doi * "}" * }
  if$
}

FUNCTION {format.pid}
{ doi empty$
{eprint empty$ 
     { url empty$ 
       { "" }
       { format.url }
     if$}
     { format.eprint }
    if$}
  { format.doi }
  if$
}

FUNCTION {n.dashify}
{
  't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
        { t #1 #2 substring$ "--" = not
            { "--" *
              t #2 global.max$ substring$ 't :=
            }
            {   { t #1 #1 substring$ "-" = }
                { "-" *
                  t #2 global.max$ substring$ 't :=
                }
              while$
            }
          if$
        }
        { t #1 #1 substring$ *
          t #2 global.max$ substring$ 't :=
        }
      if$
    }
  while$
}

FUNCTION {word.in}
{ bbl.in
  " " * }

FUNCTION {format.date}
{ year duplicate$ empty$
    { "empty year in " cite$ * "; set to ????" * warning$
       pop$ "????" }
    'skip$
  if$
  extra.label *
  before.all 'output.state :=
  after.sentence 'output.state :=
}

FUNCTION {format.btitle}
{ title
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { bbl.volume volume tie.or.space.connect
      series empty$
        'skip$
        { bbl.of space.word * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
        { series field.or.null }
        { output.state mid.sentence =
            { bbl.number }
            { bbl.number capitalize }
          if$
          number tie.or.space.connect
          series empty$
            { "there's a number but no series in " cite$ * warning$ }
            { bbl.in space.word * series * }
          if$
        }
      if$
    }
    { "" }
  if$
}

FUNCTION {is.num}
{ chr.to.int$
  duplicate$ "0" chr.to.int$ < not
  swap$ "9" chr.to.int$ > not and
}

FUNCTION {extract.num}
{ duplicate$ 't :=
  "" 's :=
  { t empty$ not }
  { t #1 #1 substring$
    t #2 global.max$ substring$ 't :=
    duplicate$ is.num
      { s swap$ * 's := }
      { pop$ "" 't := }
    if$
  }
  while$
  s empty$
    'skip$
    { pop$ s }
  if$
}

FUNCTION {convert.edition}
{ edition extract.num "l" change.case$ 's :=
  s "first" = s "1" = or
    { bbl.first 't := }
    { s "second" = s "2" = or
        { bbl.second 't := }
        { s "third" = s "3" = or
            { bbl.third 't := }
            { s "fourth" = s "4" = or
                { bbl.fourth 't := }
                { s "fifth" = s "5" = or
                    { bbl.fifth 't := }
                    { s #1 #1 substring$ is.num
                        { s eng.ord 't := }
                        { edition 't := }
                      if$
                    }
                  if$
                }
              if$
            }
          if$
        }
      if$
    }
  if$
  t
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
        { convert.edition "l" change.case$ " " * bbl.edition * }
        { convert.edition "t" change.case$ " " * bbl.edition * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
        { #1 'multiresult := }
        { t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
%        { bbl.pages pages n.dashify tie.or.space.connect }
%        { bbl.page pages tie.or.space.connect }
        { pages n.dashify }
        { pages }
      if$
    }
  if$
}

FUNCTION {first.page}
{ 't :=
  ""
    {  t empty$ not t #1 #1 substring$ "-" = not and }
    { t #1 #1 substring$ *
      t #2 global.max$ substring$ 't :=
    }
  while$
}

FUNCTION {format.journal.pages}
{ pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        {
          ", " *
          pages first.page *
        }
      if$
    }
  if$
}


FUNCTION {format.vol.num.pages}
{ volume field.or.null
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    { "" }
    { type empty$
        { bbl.chapter }
        { type "l" change.case$ }
      if$
      chapter tie.or.space.connect
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
        { word.in booktitle emphasize * }
        { word.in booktitle emphasize *
          ", " *
          editor num.names$ #1 >
            { bbl.editors }
            { bbl.editor }
          if$
          * " " *
          format.in.editors *
        }
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { bbl.techrep }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{
  word.in
  " \cite{" * crossref * "}" *
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      word.in
    }
    { bbl.volume volume tie.or.space.connect
      bbl.of space.word *
    }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{
  word.in
  " \cite{" * crossref * "}" *
}

FUNCTION {format.publisher}
{ publisher empty$
    { "empty publisher in " cite$ * warning$ }
    'skip$
  if$
  ""
  address empty$ publisher empty$ and
    'skip$
    {
      add.blank "(" *
      address empty$
        'skip$
        { address * }
      if$
      publisher empty$
        'skip$
        { address empty$
            'skip$
            { ": " * }
          if$
          publisher *
        }
      if$
      ")" *
    }
  if$
  output
}

FUNCTION {format.aaspublisher}
{ publisher empty$
    { "empty publisher in " cite$ * warning$ }
     'skip$
  if$
  ", "
  address empty$ publisher empty$ and
    'skip$
    {
      add.blank "" *
      address empty$
        'skip$
        { address * }
      if$
      publisher empty$
        'skip$
        { address empty$
            'skip$
            { ": " * }
          if$
          publisher *
        }
      if$
      "" *
    }
  if$
  output
}



STRINGS {oldname}

FUNCTION {name.or.dash}
{ 's :=
   oldname empty$
     { s 'oldname := s }
     { s oldname =
         { "---" }
         { s 'oldname := s }
       if$
     }
   if$
}

%%%%%%%%  Functions added from astrobib

FUNCTION {format.edn.btitle} % Title should be on stack.
{ duplicate$ empty$ edition empty$ or
    'skip$
    { ", " * format.edition * }
  if$
}

FUNCTION {format.ed.booktitle} % The title should be on the stack.
{ duplicate$ empty$
  { "no book title in " cite$ * warning$ "" pop$ }
  { editor empty$
    author empty$ or  % Empty author means editor already given.
      'format.edn.btitle
      { format.edn.btitle ", " * bbl.editor * " " * format.in.editors * }
      if$
  }
  if$
}

FUNCTION {format.full.book.spec} % The title should be on the stack.
{ series empty$
    { format.ed.booktitle
      volume empty$
	{ number empty$
	    'skip$
	    { " there's a number but no series in " cite$ * warning$
	      " No." number tie.or.space.connect * }
	    if$
	}
	{ ", Vol." volume tie.or.space.connect *
	  number empty$
	    'skip$
	    {"Both volume and number fields in " * cite$ * warning$ }
	    if$
	}
	if$
    }
    { volume empty$
	{ format.ed.booktitle ", " * series *
	  number empty$
	    'skip$
	    { " No." number tie.or.space.connect * }
	    if$
	}
	{ series ", Vol." volume tie.or.space.connect *  
          ", " * swap$ format.ed.booktitle *
	  number empty$
	    'skip$
	    {"Both volume and number fields in " * cite$ * warning$ }
	    if$
	}
        if$
    }
    if$
}

%%%%%%%  End of functions from astrobib

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
  crossref missing$
    { journal
      "journal" output.check
       format.vol.num.pages output
     }
     { format.article.crossref output.nonnull
       format.pages output
     }
  if$
  pages empty$
    'skip$
    { format.journal.pages }
  if$
  format.pid output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
      name.or.dash
    }
    { format.authors output.nonnull
      name.or.dash
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.date "year" output.check
  date.block
  title format.full.book.spec output
%  format.btitle "title" output.check
%  crossref missing$
%    { format.bvolume output
%      format.number.series output
%      format.edition output
   format.publisher
%    }
%    {
%      format.book.crossref output.nonnull
%    }
%  if$
  format.pages output
  format.pid output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
  format.title "title" output.check
  howpublished output
  address output
  format.pid output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
      name.or.dash
    }
    { format.authors output.nonnull
      name.or.dash
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  format.date "year" output.check
  date.block
  title format.full.book.spec output
%  format.btitle "title" output.check
%  crossref missing$
%    {
%      format.bvolume output
%      format.chapter.pages "chapter and pages" output.check
%      format.number.series output
%      format.edition output
   format.publisher
%    }
%    {
%      format.chapter.pages "chapter and pages" output.check
%      format.book.crossref output.nonnull
%    }
%  if$
  format.pages "pages" output.check
  format.pid output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
  bbl.in " " * booktitle format.full.book.spec * output
%  crossref missing$
%    { format.in.ed.booktitle "booktitle" output.check
%      format.bvolume output
%      format.number.series output
%      format.edition output
%      format.chapter.pages output
  format.publisher
%    }
%    { format.incoll.inproc.crossref output.nonnull
%      format.chapter.pages output
%    }
%  if$
  format.pages "pages" output.check
  format.pid output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output            % added
  format.date "year" output.check
  date.block
  bbl.in " " * booktitle format.full.book.spec * output
%  crossref missing$
%    { format.in.ed.booktitle "booktitle" output.check
%      format.bvolume output
%      format.number.series output
  publisher empty$
    { organization output
      address output
    }
    { organization output
      format.publisher
    }
  if$
%    }
%    { format.incoll.inproc.crossref output.nonnull
%    }
%  if$
  format.pages output
  format.pid output
  fin.entry
}

FUNCTION {conference} { inproceedings }


FUNCTION {manual}
{ output.bibitem
  format.authors output
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
  format.btitle "title" output.check
  format.edition output
  organization output
  address output
  format.pid output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
%  format.title "title" output.check
  bbl.mthesis format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.pid output
  fin.entry
}

FUNCTION {misc}
{output.bibitem
  format.authors output
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
%% new
  format.title output
%% new
  version empty$
        'skip$
 {bbl.version version * output}
 if$
%%
  howpublished output
   publisher empty$
   {}
   {format.aaspublisher}
   if$
  format.pid output
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
%  format.title "title" output.check
  bbl.phdthesis format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.pid output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization output
      organization format.key output }
    { format.editors output }
  if$
%  format.editors output
%  editor format.key output
  name.or.dash
  format.date "year" output.check
  date.block
  title format.full.book.spec output
%  format.btitle "title" output.check
%  format.bvolume output
%  format.number.series output
  publisher empty$ not      % No need for warning if no pub.
    { format.publisher }
    { editor empty$         % For empty editor, organization was already given.
        'skip$
        { organization output }
      if$
      address output
    }
  if$
%  address output
%  organization output
%  publisher output
  format.pages output
  format.pid output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
  format.title "title" output.check
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  format.pid output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output
  name.or.dash
  format.date "year" output.check
  date.block
%  format.title "title" output.check
  format.pid output
  fin.entry
}

FUNCTION {default.type} { misc }

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

FUNCTION {format.lab.names}
{ 's :=
  s #1 "{vv~}{ll}" format.name$
  s num.names$ duplicate$
  #2 >
    { pop$
      " {et~al.}" *
    }
    { #2 <
        'skip$
        { s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            {
              " {et~al.}" *
            }
            { " \& " * s #2 "{vv~}{ll}" format.name$
              * }
          if$
        }
      if$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            'key
          if$
        }
        { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {calc.short.authors}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.label
        'author.key.label
      if$
    }
  if$
  'short.list :=
}

FUNCTION {calc.label}
{ calc.short.authors
  short.list
  "("
  *
  year duplicate$ empty$
     { pop$ "????" }
     'skip$
  if$
  *
  'label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames #8 >
    { s nameptr
      "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}"
      format.name$ 't :=
      t sortify *
      "   "  *
      "zzzzz" *
    }
    { numnames 'namesleft :=
      { namesleft #0 > }
      { s nameptr
        "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}"
        format.name$ 't :=
        nameptr #1 >
          {
            "   "  *
            namesleft #1 = t "others" = and
              { "zzzzz" * }
              { t sortify * }
            if$
          }
          { t sortify * }
        if$
        nameptr #1 + 'nameptr :=
        namesleft #1 - 'namesleft :=
      }
    while$
    }
  if$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
        { "to sort, need author or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
        { key empty$
            { "to sort, need author, editor, or key in " cite$ * warning$
              ""
            }
            { key sortify }
          if$
        }
        { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.sort}
{ editor empty$
    { key empty$
        { "to sort, need editor or key in " cite$ * warning$
          ""
        }
        { key sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {presort}
{ calc.label
  label sortify
  "    "
  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
        'editor.sort
        'author.sort
      if$
    }
  if$
  #1 entry.max$ substring$
  'sort.label :=
  sort.label
  *
  "    "
  *
%% August 2016
%% commenting out this combination gets the cites right but references wrong
%%  title field.or.null
%%  sort.format.title
%%  *
%%
%% with only this line on, we get cites right but references wrong
 sort.format.title
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT

STRINGS { last.label next.extra }

INTEGERS { last.extra.num number.label }

FUNCTION {initialize.extra.label.stuff}
{ #0 int.to.chr$ 'last.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
  #0 'number.label :=
}

FUNCTION {forward.pass}
{ last.label label =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      label 'last.label :=
    }
  if$
  number.label #1 + 'number.label :=
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  extra.label 'next.extra :=
  extra.label
  duplicate$ empty$
    'skip$
    { "{\natexlab{" swap$ * "}}" * }
  if$
  'extra.label :=
  label extra.label * 'label :=
}

EXECUTE {initialize.extra.label.stuff}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {bib.sort.order}
{ sort.label
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
%% Commented out August 2016, so that multiple indentical author/year
%% combinations are entered in the order that \cite{}s are found in text.
%  title field.or.null
%  sort.format.title
%  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {bib.sort.order}

SORT

FUNCTION {begin.bib}
{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}"
  write$ newline$
  "\expandafter\ifx\csname natexlab\endcsname\relax\def\natexlab#1{#1}\fi"
  write$ newline$
  "\providecommand{\url}[1]{\href{#1}{#1}}"
  write$ newline$
  "\providecommand{\dodoi}[1]{doi:~\href{http://doi.org/#1}{\nolinkurl{#1}}}"
  write$ newline$
  "\providecommand{\doeprint}[1]{\href{http://ascl.net/#1}{\nolinkurl{http://ascl.net/#1}}}"
  write$ newline$
  "\providecommand{\doarXiv}[1]{\href{https://arxiv.org/abs/#1}{\nolinkurl{https://arxiv.org/abs/#1}}}"
  write$ newline$
}


EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
%% End of customized bst file
%%
%% End of file `aasjournal.bst'.