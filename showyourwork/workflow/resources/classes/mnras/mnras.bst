%%+
%%  Name:
%%     mnras.bst

%%  Purpose:
%%     BibTeX Bibliography style file for MNRAS

%%  Invocation:
%%     latex <latexfile>
%%     bibtex <latexfile>
%%     latex <latexfile>
%%     latex <latexfile>

%%  Description:
%%     Provides a bibliography style that will allow bibtex to be
%%     used when writing papers for MNRAS. This allows references
%%     to be inserted automatically.

%%  Notes:
%%     This file should be used in conjunction with the Monthly Notices
%%     LaTeX2e class file (mnras.cls).
%%
%%  Hyperlinks:
%%     To make the eprint and doi bibtex fields hyperlinks, you need
%%     to have \usepackage{hyperref} in your latex preamble.

%%  Authors:
%%     John P. Sleath (MRAO)
%%     Tim Jenness (MRAO/JAC)
%%     Norman Gray (Glasgow University) <norman@astro.gla.ac.uk>
%%     Oren Patashnik
%%     Susan King
%%     Sake J. Hogeveen
%%     Philip Hall
%%     Simeon Bird

%%  Bugs:
%%     The formatting has been compared to the MN style, and is
%%     believed to be compatible with it.  However this comparison has
%%     (for obvious reasons) been done most carefully for articles, so
%%     there may be deviations from it for other entry types.

%%  History:
%%     1988 (Oren Patashnik and Susan King)
%%        Apalike.bst derived from the "alpha.bst" style.
%%     1990 (Sake J. Hogeveen)
%%        Astron.bst derived from apalike.bst
%%     Summer 1995 (JPS/TJ):
%%        Extensively modified from astron.bst
%%     19 Jul 1995 (JPS/TJ):
%%        First release (version 1.1b)
%%     02 May 2001 (NG):
%%        Limit the author lists to eight people.
%%     19 Jun 2001 (TJ):
%%        Fix inproceedings
%%     30 Oct 2001 (NG):
%%        Adjust author-limiting to conform to MN style (version 1.2)
%%     22 Sep 2012 (Philip Hall):
%%        Remove comma before et al when many authors.
%%     01 July 2013 (Simeon Bird):
%%        Ported Antony Lewis' hyperlinked arxiv entry additions across.
%%        Originally modified for hyperlinked arxiv eprints by Antony Lewis April 2006
%%     August 2013 (NG):
%%        Fixed a/b/c... year suffixes -- these now appear in the
%%        correct order.
%%     February 2014 (TJ):
%%        Add MACROs for MN journal abbreviations
%%     May 2014 (NG):
%%        Shorten some functions, so that bibtex's 'single function
%%        space' parameter doesn't have to be above 100.
%%     Sep 2014 (Jared M. Gabor):
%%         Enabled ArXiv e-print identifier to be printed for articles that
%%         are pre-prints only (and have "ArXiv e-print" in the journal field).
%%         Still some problems with the use of commas.  Only 1 comma will appear
%%         between the year and "arXiv..." for such preprints, but for e.g.
%%         conference precedings there may be no comma at all.
%%     09 March 2015 (TJ):
%%         At the request of Keith Smith this file has been renamed to be mnras.bst
%%         to match the new naming convention for the MNRAS class file.
%%
%%  License:
%%     This work may be distributed and/or modified under the
%%     conditions of the LaTeX Project Public License, either version 1.3
%%     of this license or (at your option) any later version.
%%     The latest version of this license is in
%%       http://www.latex-project.org/lppl.txt
%%     and version 1.3 or later is part of all distributions of LaTeX
%%     version 2005/12/01 or later.
%%
%%     This work has the LPPL maintenance status `maintained'.
%%
%%     The Current Maintainers of this work are T. Jenness
%%     and Norman Gray <norman@astro.gla.ac.uk>
%%
%%     This work consists of the file mnras.bst

%%  This file is maintained on Github: https://github.com/timj/mn2e-bst/
%%-

ENTRY
  { address
    doi
    eprint
    archiveprefix               % used by ADS
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
%    month              not used in apalike
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
    adsurl % for MN style
    url
  }
  {}
  { label extra.label sort.label presort.key }


INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

% If threeauthorvar is set to #1 in initialise.mn2e, then the bibitem for a
% three-author reference is generated using a \mniiiauthor command.
% This command should be defined in the .cls file as follows:
%
%    \def\mniiiauthor#1#2#3{%
%      \@ifundefined{mniiiauth@#1}
%        {\global\expandafter\let\csname mniiiauth@#1\endcsname\null #2}
%        {#3}}
%
% The effect of this is that the first time that \mniiiauthor{X}{a}{b}
% is evaluated, it expands to 'a', and it thereafter expands to 'b'.
% This implements the MN style guide stipulation that "For three
% authors, give all three names at first mention, e.g. (Brown, Jones &
% Smith 1994), but use first author et al. (in roman, not italic)
% thereafter, e.g. (Brown et al. 1994)."
%
% The above definition consumes a command for each such author.  It's
% surely possible to avoid this (with some constructions involving
% {\\#1}; see Appendix D cleverness), but that would verge on the
% arcane, and not be really worth it.
INTEGERS { maxauthors threeauthorvar }

FUNCTION {initialise.mn2e}
{
  #8 'maxauthors :=
  #0 'threeauthorvar :=
}

STRINGS { s t fmt sep }

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { ". " *  write$
          newline$
%          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { ", " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

% output.nonnull.extra is like output.nonnull, except that the thing
% we output is 'extra' to a block, so if we're in the after.block
% state, we append it rather than starting a new block.
FUNCTION {output.nonnull.extra}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
        { ", " *  write$
          newline$
          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { ", " * write$ }
          if$
        }
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output.nonnull.nocomma}
{'s :=
  output.state mid.sentence =
    { write$ }
    { output.state after.block =
        { ". " *  write$
          newline$
%          "\newblock " write$
        }
        { output.state before.all =
            'write$
            { ", " * write$ }
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

FUNCTION {output.extra}
{ duplicate$ empty$
    'pop$
    'output.nonnull.extra
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

FUNCTION {output.check.nocomma}
{ 't :=
  duplicate$ empty$
    'pop$ % or warn with { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull.nocomma
  if$
}

FUNCTION {output.check.extra}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull.extra
  if$
}

FUNCTION {output.year.check}
{ year empty$
    { "empty year in " cite$ * warning$ }
    { write$
      " " year * extra.label *
      mid.sentence 'output.state :=
    }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem[" write$
  label write$
  "]{" write$
  cite$ write$
  "}" write$
  % debugging: define eg \def\logsortkey#1{{[\tiny #1]}} to see this
  %"\logsortkey{" presort.key * "/" * sort.key$ * "}" * write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ write$
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

FUNCTION {boldface}
{ duplicate$ empty$
    { pop$ "" }
    { "{ " swap$ * "}" * }
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{ " swap$ * "}" * }
  if$
}

% Make hyperref links.
% The LaTeX hyperref package supports a macro \href{URL}{text}
%
% On call, the stack contains
%   href-url
%   link-text
FUNCTION {make.href.hyperref}
{
  "\href {" swap$ * "} {" * swap$ * "}" *
}

INTEGERS { nameptr namesleft numnames }

% Helper functions for print.all.names: print the last or 'others' name.
%
% On call, the stack contains
%  a formatted name,
%  and the built-up string of names so far.
FUNCTION {print.all.names.plain.trailer}
{ duplicate$ "others" =
    { pop$ " et~al.," * }
    { " " swap$ * * }
  if$
}

% as above
FUNCTION {print.all.names.label.trailer}
{ duplicate$ "others" =
    { pop$ " et~al." * }
    { " \& " swap$ * * }
  if$
}

% Print out the complete author list
% (because we have maxauthors or fewer).
%
% On call, the stack contains
%   a format string (for format.name$),
%   a name-separator string (such as ", "),
%   an author or editor field,
%   and one of the 'trailer' functions above.
FUNCTION {print.all.names}
{ 'fmt :=
  'sep :=
  's :=
  s num.names$ 'namesleft :=
  namesleft #2 < % we're not going to need the trailer function, so discard it
    'pop$
    'skip$
  if$
  #1 'nameptr :=
  { namesleft #0 > }
    { s nameptr fmt format.name$ 't := % last name first
      nameptr #1 >
        { namesleft #1 >
            { sep * t * }
            { % namesleft == 1: this is the last name.
              % Invoke the trailer function.
              s num.names$ #2 > % _why_ do we add this extra space for long lists?
                { " " * }
                'skip$
              if$
              % This next line is a bit clever.
              %
              % We start with a stack
              %   string containing the list of names so far
              %   trailer function
              %
              % Rearrange the stack so that we have
              %   'skip
              %   trailer function
              %   #1
              %   t
              %   names so far
              % then call if$.  This is a rather indirect way of
              % invoking the trailer function on the stack
              %   t
              %   names so far
              %
              % Ta dahhhh!
              swap$ t swap$ #1 swap$ 'skip$ if$
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

% Name formatting.  According to the MN instructions for authors, `In
% the references list [...] List all authors if eight or fewer,
% otherwise et al.'.  This means (according to a ruling from the MN
% editors) that if there are eight or more authors, they should be
% listed as "first one plus et al."
%
% On call, the stack contains
%   an author or editor field.
FUNCTION {format.names}
{ duplicate$ num.names$ maxauthors >
    { % Too many authors -- format as First author et al.
      #1 "{vv }{ll}{ jj}{ f., }" format.name$  % last name first
      " et~al.," *
    }
    { 'print.all.names.plain.trailer swap$ " " "{vv }{ll}{ jj}{ f., }" print.all.names }
  if$
}

% Format an author/editor field, for inclusion in a label
% (ie, "Surname, Surname \& Surname")
%
% On call, the stack contains an author/editor field.
FUNCTION {format.names.label}
{ duplicate$ num.names$ maxauthors >
    { % Too many authors -- format as First author et al.
      #1 "{vv }{ll}{ jj}" format.name$  % last name first
      " et~al.," *
    }
    { 'print.all.names.label.trailer swap$ ", " "{vv }{ll}" print.all.names }
  if$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.key}                   % this function is just for apalike
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
        { " eds," * }
        { " ed." * }
      if$
    }
  if$
}

FUNCTION {format.editors.extra}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
        { " eds" * }
        { " ed." * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }
    { title "t" change.case$ }
  if$
}

FUNCTION {format.eprint}
{ eprint empty$
    { "" }
    { " (\mn@eprint {"
      archiveprefix empty$
        { "" }
        { archiveprefix }
      if$
      * "} {" * eprint * "})" * }
    if$
}

FUNCTION {format.doi}
{ doi empty$
    { "" }
    { ", \mn@doi{" * doi * "}" }
    if$
}

FUNCTION {n.dashify}
{ 't :=
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

FUNCTION {first.page.number}
{ 't :=
  ""
    { t "" =
        { #0 }
        { t #1 #1 substring$ "-" = not }
      if$
    }
    { t #1 #1 substring$ *
      t #2 global.max$ substring$ 't :=
    }
  while$
}

FUNCTION {format.btitle}
{ title
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #5 <
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
    { "" series tie.or.space.connect
      volume empty$
        'skip$
        { " Vol. " * volume * }
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
            { "No." }
            { "No." }
          if$
          number tie.or.space.connect
          series empty$
            { "there's a number but no series in " cite$ * warning$ }
            { " in " * series * }
          if$
        }
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
        { edition "l" change.case$ " edn" * }
        { edition "t" change.case$ " edn" * }
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
        { "pp" pages n.dashify tie.or.space.connect }
        { "p." pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.page}
{ pages empty$
    { "" }
    { "p.~" pages first.page.number * }
  if$
}

% If a paper is available as a preprint only, ADS will call the
% journal "ArXiv e-print".  In this case MN suggests the journal
% should be called "preprint".
FUNCTION {format.journal.doi}
{ doi empty$
    { journal "ArXiv e-prints" =
        { "preprint" }
        { journal }
      if$
    }
    { journal "ArXiv e-prints" =
        { "preprint" }
        { "\mn@doi [" journal * "] {" * doi * "}" * }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
{ volume field.or.null
  number empty$
    'skip$
    { "(" number * ")" * *
      volume empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { ", " * pages n.dashify * }
      if$
    }
  if$
}

FUNCTION {format.vol.page}
{ volume field.or.null
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { ", " * pages first.page.number * }
      if$
    }
  if$
}

FUNCTION {format.vol.num.page}
{ volume field.or.null
  number empty$
    'skip$
    { "(" number * ")" * *
      volume empty$
        { "there's a number but no volume in " cite$ * warning$ }
        'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
        { pop$ format.pages }
        { ", " * pages first.page.number * }
      if$
    }
  if$
}
FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
        { "Chapt." }
        { type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
        'skip$
        { ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
        { "in " booktitle * }
	{ format.bvolume empty$
      		{ "in " format.editors.extra * ", " * booktitle * }
        	{ "in " format.editors.extra * ", " * format.bvolume * ", " * booktitle * }
		if$
	}
      if$
    }
  if$
}

FUNCTION {format.in.ed.collect}
{ booktitle empty$
    { "" }
    { editor empty$
      { "in " }
      { "in " format.editors.extra * ", " * }
      if$
      series empty$
      { skip$ }
      { series * }
      if$
      volume empty$
      { number empty$
        { }
        { output.state mid.sentence =
            { "No." }
            { "No." }
          if$
          number tie.or.space.connect *
        }
        if$
      }
      {
	output.state mid.sentence =
        { ", " * }
	{ skip$ }
	if$
        "Vol." volume tie.or.space.connect *
      }
      if$
      output.state mid.sentence =
      { ", " * }
      { skip$ }
      if$
      booktitle *
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

% the \url{} macro is supplied by hyperref
FUNCTION {format.url}
{
  url empty$
    { "" }
    { "\url {" url * "}" * }
  if$
}

FUNCTION {format.tr.number}
{
  type empty$
    { "Technical Report" }
    { type }
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ "in"                                                  % this is for apalike
  " \cite{" * crossref * "}" *
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "in "
    }
    { "Vol." volume tie.or.space.connect
      " of " *
    }
  if$
  "\cite{" * crossref * "}" *                           % this is for apalike
}

FUNCTION {format.incoll.inproc.crossref}
{ "in"                                                  % this is for apalike
  " \cite{" * crossref * "}" *
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  % MN style does not display eprints or DOIs on published articles,
  % but does link the journal name to any available DOI.
  crossref missing$
    { format.journal.doi "journal" output.check.extra
      % ...and links volume and page information to any available adsurl
      adsurl empty$
        { format.vol.page }
        { format.vol.page adsurl make.href.hyperref }
      if$
      output.extra
    }
    { format.article.crossref output.nonnull
      format.page output
    }
  if$

  % If the "journal" field is the arXiv, then print out the arXiv identifier
  journal "ArXiv e-prints" =
    { format.eprint "eprint" output.check.nocomma }
    { }
  if$

  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  output.year.check                             % special for apalike
  format.btitle "title" output.check
  format.edition output
  new.block
  crossref missing$
    { format.bvolume output
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
    }
      {format.book.crossref output.nonnull}
  if$
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output                              % special for
  output.year.check                                     % apalike
  format.btitle "title" output.check
  howpublished output
  address output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check
      editor format.key output
    }
    { format.authors output.nonnull
      crossref missing$
        { "author and editor" editor either.or.check }
        'skip$
      if$
    }
  if$
  output.year.check                             % special for apalike
  format.btitle "title" output.check
  format.edition output
  new.block
  crossref missing$
    {
      new.sentence
      publisher "publisher" output.check
      address output
      format.pages output
    }
    { format.chapter.pages "chapter and pages" output.check
      format.book.crossref output.nonnull
    }
  if$
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  crossref missing$
    {
      format.in.ed.collect "booktitle" output.check
      format.edition output
      booktitle empty$ edition empty$ =
        { } { new.block } if$
      new.sentence
      publisher "publisher" output.check
      address output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  format.chapter.pages output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  crossref missing$
    {

      format.in.ed.booktitle "booktitle" output.check
      ". " *
%      format.bvolume output
      before.all 'output.state :=
      new.block
%      format.title "title" output.check
      format.number.series output
      new.block
%      organization output
      publisher output
      address output
      format.pages output
    }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  format.authors output
  author format.key output                              % special for
  output.year.check                                     % apalike
  format.btitle "title" output.check
  new.block
  organization address new.block.checkb
  organization output
  address output
  format.edition output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  author format.key output                              % special for apalike
  output.year.check
  format.btitle output
  howpublished output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  format.editors output
  editor format.key output                              % special for
      mid.sentence 'output.state :=
  output.year.check                                     % apalike
  format.btitle "title" output.check
  " " *
  before.all 'output.state :=
%  new.block
  format.bvolume output
  format.number.series output
  new.block
%  organization output                           % a nonempty organization
  publisher output
  address output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  format.tr.number output.nonnull
  format.btitle "title" output.check
  format.url output
  new.block
  institution "institution" output.check
  address output
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for
  output.year.check                                     % apalike
  format.btitle "title" output.check
  note "note" output.check
  % the following are presumably redundant for 'unpublished'
  format.eprint "eprint" output.check.nocomma
  format.doi "doi" output.check.nocomma
  format.url output
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sept."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

%% Journal macros for bib files. Should be kept synchronized with
%% the entries in the class file.

%% Astronomy and Astrophysics
MACRO {aap} {"A\&A"}

%% Astronomy and Astrophysics
MACRO {astap} {"A\&A"}

%% Astronomy and Astrophysics Review (the)
MACRO {aapr} {"A\&ARv"}

%% Astronomy and Astrophysics Supplement Series
MACRO {aaps} {"A\&AS"}

%% Acta Astronomica
MACRO {actaa} {"Acta Astron."}

%% Astrofizika
MACRO {afz} {"Afz"}

%% Astronomical Journal (the)
MACRO {aj} {"AJ"}

%% Applied Optics
MACRO {ao} {"Appl. Opt."}

%% Applied Optics
MACRO {applopt} {"Appl. Opt."}

%% Astrophysics Letters
MACRO {aplett} {"Astrophys.~Lett."}

%% Astrophysical Journal
MACRO {apj} {"ApJ"}

%% Astrophysical Journal, Letters
MACRO {apjl} {"ApJ"}

%% Astrophysical Journal, Letters
MACRO {apjlett} {"ApJ"}

%% Astrophysical Journal, Supplement
MACRO {apjs} {"ApJS"}

%% Astrophysical Journal, Supplement
MACRO {apjsupp} {"ApJS"}

%% Astrophysics and Space Science
MACRO {apss} {"Ap\&SS"}

%% Annual Review of Astronomy and Astrophysics
MACRO {araa} {"ARA\&A"}

%% Astronomicheskii Zhurnal
MACRO {azh} {"Azh"}

%% Bulletin of the American Astronomical Society
MACRO {baas} {"BAAS"}

%% Bulletin of the Astronomical Institutes of Czechoslovakia
MACRO {bac} {"Bull. Astron. Inst. Czechoslovakia"}

%% Bulletin Astronomical Institute of the Netherlands
MACRO {bain} {"Bull. Astron. Inst. Netherlands"}

%% Chinese Astronomy and Astrophysics
MACRO {caa} {"Chinese Astron. Astrophys."}

%% Chinese Journal of Astronomy and Astrophysics
MACRO {cjaa} {"Chinese J.~Astron. Astrophys."}

%% Fundamentals of Cosmic Physics
MACRO {fcp} {"Fundamentals Cosmic Phys."}

%% Geochimica Cosmochimica Acta
MACRO {gca} {"Geochimica Cosmochimica Acta"}

%% Geophysics Research Letters
MACRO {grl} {"Geophys. Res. Lett."}

%% IAU Cirulars
MACRO {iaucirc} {"IAU~Circ."}

%% Icarus
MACRO {icarus} {"Icarus"}

%% Journal of Astrophysics and Astronomy
MACRO {japa} {"J.~Astrophys. Astron."}

%% Journal of Cosmology and Astroparticle Physics
MACRO {jcap} {"J.~Cosmology Astropart. Phys."}

%% Journal of Chemical Physics
MACRO {jcp} {"J.~Chem.~Phys."}

%% Journal of Geophysics Research
MACRO {jgr} {"J.~Geophys.~Res."}

%% Journal of Quantitiative Spectroscopy and Radiative Transfer
MACRO {jqsrt} {"J.~Quant. Spectrosc. Radiative Transfer"}

%% Journal of the RAS of Canada
MACRO {jrasc} {"J.~R.~Astron. Soc. Canada"}

%% Memoirs of the RAS
MACRO {memras} {"Mem.~RAS"}

%% Memoire della Societa Astronomica Italiana
MACRO {memsai} {"Mem. Soc. Astron. Italiana"}

%% Monthly Notes of the Astronomical Society of Southern Africa
MACRO {mnassa} {"MNASSA"}

%% Monthly Notices of the Royal Astronomical Society
MACRO {mnras} {"MNRAS"}

%% New Astronomy
MACRO {na} {"New~Astron."}

%% New Astronomy Review
MACRO {nar} {"New~Astron.~Rev."}

%% Nature
MACRO {nat} {"Nature"}

%% Nuclear Physics A
MACRO {nphysa} {"Nuclear Phys.~A"}

%% Physical Review A: General Physics
MACRO {pra} {"Phys. Rev.~A"}

%% Physical Review B: Solid State
MACRO {prb} {"Phys. Rev.~B"}

%% Physical Review C
MACRO {prc} {"Phys. Rev.~C"}

%% Physical Review D
MACRO {prd} {"Phys. Rev.~D"}

%% Physical Review E
MACRO {pre} {"Phys. Rev.~E"}

%% Physical Review Letters
MACRO {prl} {"Phys. Rev.~Lett."}

%% Publications of the Astronomical Society of Australia
MACRO {pasa} {"PASA"}

%% Publications of the Astronomical Society of the Pacific
MACRO {pasp} {"PASP"}

%% Publications of the Astronomical Society of Japan
MACRO {pasj} {"PASJ"}

%% Physics Reports
MACRO {physrep} {"Phys.~Rep."}

%% Physica Scripta
MACRO {physscr} {"Phys.~Scr."}

%% Planetary Space Science
MACRO {planss} {"Planet. Space~Sci."}

%% Proceedings of the Society of Photo-Optical Instrumentation Engineers
MACRO {procspie} {"Proc.~SPIE"}

%% Revista Mexicana de Astronomia y Astrofisica
MACRO {rmxaa} {"Rev. Mex. Astron. Astrofis."}

%% Quarterly Journal of the RAS
MACRO {qjras} {"QJRAS"}

%% Science
MACRO {sci} {"Science"}

%% Sky and Telescope
MACRO {skytel} {"Sky \& Telesc."}

%% Solar Physics
MACRO {solphys} {"Sol.~Phys."}

%% Soviet Astronomy (aka Astronomy Reports)
MACRO {sovast} {"Soviet~Ast."}

%% Space Science Reviews
MACRO {ssr} {"Space Sci. Rev."}

%% Zeitschrift fuer Astrophysik
MACRO {zap} {"Z.~Astrophys."}

%% Astronomical Society of the Pacific Conference Series
MACRO {aspconf} {"ASP Conf. Ser."}

%% Astronomische Nachrichten
MACRO {astnach} {"Astron. Nachr."}

%% Astronomy & Computing
MACRO {astcomp} {"Astron. Comput."}

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

% There are three apalike cases: one person (Jones),
% two (Jones and de~Bruijn), and more (Jones et~al.).
% This function is much like format.crossref.editors.
%
% The goal here is to be compatible with the MN style prescription that:
%     If there are several references with the same first author,
%     arrange in the following order:
%     firstly single-author papers (by date);
%     then two-author papers (alphabetically by co-author, then by date);
%     then multi-author papers (by date).
% We deem 'Jones and others' to be a multi-author paper, in this sense,
% so that it sorts, by date, with the three-or-more author papers.
% An alternative is to deem 'Jones and others' to sort just after the single-author papers
% (this case isn't really covered by the above rules, but this position looks sensible);
% we can get that by putting extra spaces before the 'et~al.', as noted below.
% Note: we can't put arbitrary characters here, to force the order,
% because the result of this function appears as the citation label in the typeset result.
FUNCTION {format.lab.names}
{ 's :=
  s #1 "{vv~}{ll}" format.name$
  s num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
        'skip$
        { s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
            { " et~al." * } % <-----+-- alternate sortings
            %{ "   et~al." * } % <--+
            { " \& " * s #2 "{vv~}{ll}" format.name$ * }
          if$
        }
      if$
    }
  if$
}

% Format names in full
%
FUNCTION {format.lab.names.full}
{ 's :=
  s #1 "{vv~}{ll}" format.name$
  s num.names$ duplicate$
  #2 <
    'skip$
    { s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	'skip$
        { " \& " * s #2 "{vv~}{ll}" format.name$ * }
      if$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key                                    % apalike uses the whole key
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
            'key                                % apalike uses the whole key
          if$
        }
        { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {author.editor.key.label.full}
{ author empty$
    { editor empty$
        { key empty$
            { cite$ #1 #3 substring$ }
            'key
          if$
        }
        { editor format.names.label }
      if$
    }
    { author format.names.label }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
        { cite$ #1 #3 substring$ }
        'key                    % apalike uses the whole key, no organization
      if$
    }
    { editor format.lab.names }
  if$
}

% Calculate the 'label' for the entry, which is used for sorting.
% Note that this does not include the a/b/c... for disambiguating
% years (that's added later, in forward.pass).
FUNCTION {calc.label}
{ % Define t: this is the 'full' set of names
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label.full
    { type$ "proceedings" =
        { editor format.names.label }
        { author format.names.label }
      if$
    }
  if$
  't :=
  % Define s: this is the 'shortened' set of names, with 'et al.',
  % which are to appear in citations
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.label                       % apalike ignores organization
        'author.key.label                       % for labeling and sorting
      if$
    }
  if$
  's :=
  % the '\protect... is for apalike
  "\protect\citeauthoryear{" t * "}{" *
  threeauthorvar
  numnames #3 =
  and
    { "\protect\mniiiauthor{" cite$ * "}{" *
      t * "}{" * s * "}" * }
    { s }
  if$
  *
  "}{" *
  year field.or.null purify$ #-1 #4 substring$  % use all four digits
  *
  % the matching closing "}" comes in at the reverse.pass
  'label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
        { "   " * }
        'skip$
      if$                                               % apalike uses initials
      s nameptr "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't := % <= here
      "Mc" t #2 text.prefix$ =
          { "mac" t #3 t text.length$ substring$ * 't := }
          'skip$
      if$
      nameptr numnames = t "others" = and
        { " et~al." * }
        { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
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

% Format numbers (eg volume, number, pages) into a fixed-width string,
% for sorting.
FUNCTION {number.for.presort}
{ field.or.null
  duplicate$ text.length$ 'len :=
  ""
  { #5 len > }
    { "0" *
      len #1 + 'len := }
  while$
  swap$ *
}

% apalike uses two sorting passes; the first one sets the
% labels so that the `a's, `b's, etc. can be computed;
% the second pass puts the references in "correct" order.
% The presort function is for the first pass. It computes
% presort.key, sort.label, and sort.key$
FUNCTION {presort}
{ calc.label
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
        'editor.key.label       % apalike ignores organization
        'author.key.label       % for labeling and sorting
      if$
    }
  if$
  purify$                       % remove {braces} round names, to avoid confusing sorting
  "    "
  *
  year number.for.presort *
  #1 entry.max$ substring$
  'presort.key :=
  presort.key
  "    "
  *
  % Style apalike prefers 'title field.or.null sort.format.title *' here.
  % We don't, because we will predominantly be sorting articles in
  % this style, and we don't include their titles in the final bibliography.
  volume number.for.presort *
  number number.for.presort *
  pages empty$
    { "0" }
    { pages first.page.number }
  if$
  number.for.presort *
  #1 entry.max$ substring$
  'sort.key$ :=
  sort.key$ 'sort.label :=
}

EXECUTE {initialise.mn2e}

ITERATE {presort}

SORT            % by label, sort.label, title---for final label calculation

% Manage the construction of the 'extra' labels for references.  These
% are the a/b/c/... year suffixes which distinguish different
% publications by the same author.
%
% These apalike labels are only for the text; the 'label' does not
% appear in the bibiliography, so that the 'longest.label' variable is
% unnecessary in this style.
STRINGS { last.key next.extra }

INTEGERS { last.extra.num }

FUNCTION {initialize.extra.label.stuff}
{ #0 int.to.chr$ 'last.key :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

% By the time that this function is iterated over, the entries have
% had presort.key assigned, and have been presorted.  However these labels
% do not (by design) respect the presence of multiple entries in the
% same year, so that there are duplicates.  Detect successive
% duplicates, and add an extra a/b/c... suffix.
%
% At the end of this forward.pass, each entry with a duplicate presort.key
% will have an extra.label, starting with 'b' for the first duplicate,
% 'c' for the second, and so on.  That is, only the first
% soon-to-be-discovered-to-be-duplicate label will not have an
% extra.label.
FUNCTION {forward.pass}
{ last.key presort.key =
    { last.extra.num #1 + 'last.extra.num :=
      last.extra.num int.to.chr$ 'extra.label :=
    }
    { % first time we've seen this key
      "a" chr.to.int$ 'last.extra.num :=
      "" 'extra.label :=
      presort.key 'last.key :=
    }
  if$
}

% Work back through the labels.  If, on this pass we discover that the
% next entry in the sequence (which is the last entry processed in
% this reverse pass) had an extra.label of 'b', then this tells us
% that _this_ entry should have an extra.label of 'a', which we were
% unable to detect in the forward.pass above.
%
% We also close the brace pair {...} opened in calc.label
FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
    'skip$
  if$
  label extra.label * "}" * 'label :=
  extra.label 'next.extra :=
}

EXECUTE {initialize.extra.label.stuff}

ITERATE {forward.pass}

REVERSE {reverse.pass}

% Now that the label is right we sort for real.  We sort on sort.label
% then year then title.  Since sort.label starts with presort.key,
% this sorting is compatible with the letter suffixes added to years.
%
% This is for the second sorting pass.
FUNCTION {bib.sort.order}
{ sort.label
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {bib.sort.order}

SORT            % by sort.label, year, title---giving final bibliography order

FUNCTION {begin.bib}
{
  preamble$ empty$                              % no \etalchar in apalike
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}" write$ newline$ % no labels in apalike
  "\makeatletter" write$ newline$
  % The following is the contents of mnras-insert.sty, minus comments and whitespace.
  % DO NOT EDIT THIS TEXT,
  % but instead make mnras-insert.txt, and paste the results in here.
"\relax" write$ newline$
"\def\mn@urlcharsother{\let\do\@makeother \do\$\do\&\do\#\do\^\do\_\do\%\do\~}" write$ newline$
"\def\mn@doi{\begingroup\mn@urlcharsother \@ifnextchar [ {\mn@doi@} {\mn@doi@[]}}" write$ newline$
"\def\mn@doi@[#1]#2{\def\@tempa{#1}\ifx\@tempa\@empty \href {http://dx.doi.org/#2} {doi:#2}\else \href {http://dx.doi.org/#2} {#1}\fi \endgroup}" write$ newline$
"\def\mn@eprint#1#2{\mn@eprint@#1:#2::\@nil}" write$ newline$
"\def\mn@eprint@arXiv#1{\href {http://arxiv.org/abs/#1} {{\tt arXiv:#1}}}" write$ newline$
"\def\mn@eprint@dblp#1{\href {http://dblp.uni-trier.de/rec/bibtex/#1.xml} {dblp:#1}}" write$ newline$
"\def\mn@eprint@#1:#2:#3:#4\@nil{\def\@tempa {#1}\def\@tempb {#2}\def\@tempc {#3}\ifx \@tempc \@empty \let \@tempc \@tempb \let \@tempb \@tempa \fi \ifx \@tempb \@empty \def\@tempb {arXiv}\fi \@ifundefined {mn@eprint@\@tempb}{\@tempb:\@tempc}{\expandafter \expandafter \csname mn@eprint@\@tempb\endcsname \expandafter{\@tempc}}}" write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\makeatother" write$ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}