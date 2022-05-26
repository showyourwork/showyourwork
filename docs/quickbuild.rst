.. admonition:: Building articles with |showyourwork|

   Find an open source article you'd like to reproduce or tinker with?
   Build it in 3 easy steps.

   1. Install the latest version of |showyourwork|:
   
      .. code-block:: bash
      
         pip install -U showyourwork

      .. raw:: html

        <br/>

   2. Clone the article repository (replace ``user`` and ``repo`` with the GitHub user name
      and repository name for the desired project):
   
      .. code-block:: bash
      
         git clone https://github.com/user/repo
         cd repo
      
      .. raw:: html

        <br/>

   3. Build!

      .. code-block:: bash
      
         showyourwork

      .. raw:: html

         <pre>
         <span style="color:green;">Setting up the workflow...</span>
         <span style="color:green;">Generating the article PDF...</span>
         <span style="color:green;">Done!</span>
         </pre>

   After installing all the required packages and running any pipeline or
   figure scripts, this will generate the article PDF ``ms.pdf`` in the current 
   working directory.