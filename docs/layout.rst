Repository layout
=================


.. raw:: html

    <style>
      /*
           https://codepen.io/asraven/pen/qbrQMX
      */
      .directory-list ul {
        margin-left: 10px;
        padding-left: 20px;
        border-left: 1px dashed #ddd;
      }
      .directory-list li {
        list-style: none;
        color: #888;
        font-size: 17px;
        font-style: normal;
        font-weight: normal;
      }
      .directory-list a {
        border-bottom: 1px solid transparent;
        color: #888;
        text-decoration: none;
        transition: all 0.2s ease;
      }
      .directory-list a:hover {
        border-color: #eee;
        color: #000;
      }
      .directory-list .folder,
      .directory-list .folder > a {
        color: #777;
        font-weight: bold;
      }
      .directory-list li:before {
        margin-right: 10px;
        content: "";
        height: 20px;
        vertical-align: middle;
        width: 20px;
        background-repeat: no-repeat;
        display: inline-block;
        /* file icon by default */
        background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><path fill='lightgrey' d='M85.714,42.857V87.5c0,1.487-0.521,2.752-1.562,3.794c-1.042,1.041-2.308,1.562-3.795,1.562H19.643 c-1.488,0-2.753-0.521-3.794-1.562c-1.042-1.042-1.562-2.307-1.562-3.794v-75c0-1.487,0.521-2.752,1.562-3.794 c1.041-1.041,2.306-1.562,3.794-1.562H50V37.5c0,1.488,0.521,2.753,1.562,3.795s2.307,1.562,3.795,1.562H85.714z M85.546,35.714 H57.143V7.311c3.05,0.558,5.505,1.767,7.366,3.627l17.41,17.411C83.78,30.209,84.989,32.665,85.546,35.714z' /></svg>");
        background-position: center 2px;
        background-size: 60% auto;
      }
      .directory-list li.folder:before {
        /* folder icon if folder class is specified */
        background-image: url("data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 100 100'><path fill='lightblue' d='M96.429,37.5v39.286c0,3.423-1.228,6.361-3.684,8.817c-2.455,2.455-5.395,3.683-8.816,3.683H16.071 c-3.423,0-6.362-1.228-8.817-3.683c-2.456-2.456-3.683-5.395-3.683-8.817V23.214c0-3.422,1.228-6.362,3.683-8.817 c2.455-2.456,5.394-3.683,8.817-3.683h17.857c3.422,0,6.362,1.228,8.817,3.683c2.455,2.455,3.683,5.395,3.683,8.817V25h37.5 c3.422,0,6.361,1.228,8.816,3.683C95.201,31.138,96.429,34.078,96.429,37.5z' /></svg>");
        background-position: center top;
        background-size: 75% auto;
      }
    </style>

      <div class="box">
        <ul class="directory-list">
          <li class="folder">.github
            <ul>
              <li class="folder">workflows
                <ul>
                  <li>showyourwork.yml</li>
                </ul>
              </li>
            </ul>
          </li>
          <li class="folder">showyourwork
          </li>
          <li class="folder">src
            <ul>
              <li class="folder">figures
                <ul>
                  <li>.gitignore</li>
                  <li>matplotlibrc</li>
                </ul>
              </li>
              <li class="folder">static
                <ul>
                  <li>.gitignore</li>
                </ul>
              </li>
              <li>.gitignore</li>
              <li>bib.bib</li>
              <li>ms.tex</li>
            </ul>
          </li>
          <li>.gitignore</li>
          <li>.gitmodules</li>
          <li>LICENSE</li>
          <li>README.md</li>
          <li>Snakefile</li>
          <li>environment.yml</li>
        </ul>
      </div>
