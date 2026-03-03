#!/usr/bin/env python3

"""
This script creates a html report by reading ReadME.txt and images from folders/subfolders of a project

Parameters
----------
basedir : str
    Path to base directory where the report data is stored
readme_file : str
    Path to txt report

Instructions
------------
The report has to be written as a plain txt file with specific sections/elements defined by headers (see below)

Headers
-------
HEADER
    Syntax:
        HEADER:<header_name>
    e.g.
        HEADER:1 - Data pre-processing

SUBHEADER
    Syntax:
        SUBHEADER:<subheader_name>
    e.g.
        SUBHEADER:1.1 - Data normalization and filtering

PLOT
    Syntax:
        PLOT:<plot_header>:<relative_path_to_plot>
    e.g.
        PLOT:Cd3e expression:./GenesOfInterest/Cd3e_Expression_UMAP-30.png

TABLE
    Syntax:
        TABLE:<table_header>:<relative_path_to_tab_separated_table>
    e.g.
        TABLE:Predicted cluster identities:./Clusters/ClusterIdentities.tsv
"""

### ---------------------------------------- ###

class produceReport:
    
    def __init__(self, readme_file):
        
        self.blocks = open(readme_file, 'r').read().split('\n\n')
        
        self.processComponents()
    
    ### ------------------------------------ ###

    def processComponents(self):
        
        # Init html file
        html = ['<!DOCTYPE html>', '<html>']
        
        # Adding navbar
        navbar = self.makeNavbar()
        html.append(navbar)
        
        # Adding folder structure
        folder_structure = self.getFolderStructure().replace('\n', '<br>').replace(' ', '&nbsp')
        
        new_element = f"""
        <h1 id="ch0" style="color: white; background-color: #4c4cff; border: 2px solid black; border-radius: 5px 5px 5px 5px;">Folder structure</h1>
        <body>
            <p style="line-height: 1.5; font-size: 20px;font-family: monospace;"><b>{folder_structure}</b></p>
        </body>
        """
        
        html.append(new_element)
        
        # Processing each component
        self.header_counts = 1
        for block in self.blocks:
            
            paragraphs = block.split('\n')
            new_section = []
            
            # Parsing body
            new_section.append(self.makeBody(paragraphs))
            
            new_section = '\n'.join(new_section)
            html.append(new_section)
        
        # Completing html file
        html.append('</html>')
        html = '\n'.join(html)
        
        # Saving html
        with open('analysis_report.html', 'w') as output:
            output.write(html)
        
        self.html = html

    ### ------------------------------------ ###

    def makeNavbar(self):
        
        headers = ['HEADER:Folder structure'] + [line for b in self.blocks for line in b.split('\n') if 'HEADER:' in line or 'SUBHEADER:' in line]

        # Init navbar
        navbar = """
        <head>
        <style>

        #navbar {
        float: none;
        overflow: scroll;
        width: 500px;
        min-height: 250px;
        max-height: 500px;
        position: fixed;
        right: -450px;
        border: 2px solid black;
        border-radius: 5px 5px 5px 5px;
        }
    
        #navbar:hover{
        right: 0px;
        background-color: white;
        }

        #navbar a {
        position: absolute;
        left: 0px;
        color: white;
        text-decoration: none;
        padding: 2px;
        border: 2px solid black;
        border-radius: 5px 5px 5px 5px;
        opacity: 0.5;
        white-space: nowrap;
        }
        
        """

        # Formatting individual navbar entries in <style>
        for n in range(len(headers)):
            
            if headers[n][:7] == 'HEADER:':
                
                navbar += f"""
                #ch{n}_nav {{
                top: {32 * n}px;
                font-size: 16px;
                background-color: #6a6a6a;
                }}
    
                #ch{n}_nav:hover {{
                opacity: 1;
                background-color: #4c4cff;
                }}
                """
            
            else:
                
                navbar += f"""
                #ch{n}_nav {{
                top: {32 * n}px;
                font-size: 16px;
                background-color: #aaa;
                }}
    
                #ch{n}_nav:hover {{
                opacity: 1;
                background-color: #4c97ff;
                }}
                """
        
        # Closing <style>
        navbar += """
        
        </style>
        </head>
        """
        
        # Creating navbar div object
        navbar += """
        <div id="navbar">
        """
        
        for n,h in enumerate(headers):
            
            if h[:7] == 'HEADER:':
                
                navbar += f"""
                <a href="#ch{n}" id="ch{n}_nav">{h.replace('HEADER:', '')}</a>
                """
            
            else:
                
                navbar += f"""
                <a href="#ch{n}" id="ch{n}_nav">{h.replace('SUBHEADER:', '')}</a>
                """
        
        navbar += """
        </div>
        """
        
        return navbar
    
    ### ------------------------------------ ###

    @staticmethod
    def getFolderStructure():
        
        # Finding directories and their level in the hierarchy
        directories = [directory[0] for directory in walk(top='./')]
        dir_levels = [current_dir.count('/') - 1 for current_dir in directories]
        dir_levels[0] = -1
        
        # Building structure
        structure = []
        n_spaces = 4
        for n, current_dir in enumerate(directories):
            
            if n == 0: # First element (root)
                
                structure.append(current_dir)
            
            else:
                
                # Extracting info
                basename = current_dir.split('/')[-1]
                current_level = dir_levels[n]
                next_element_level = dir_levels[n + 1] if n + 1 != len(directories) else dir_levels[n] - 1
                
                # Building new line of the structure
                line_elements = []
                for l in range(current_level + 1):
                    
                    if l < dir_levels[n]:
                        
                        if l in dir_levels[n:]:
                            
                            new_element = f"\u2502{' ' * (n_spaces - 1)}"
                        
                        else:
                            
                            new_element = f"{' ' * n_spaces}"
                        
                    else:
                        
                        if current_level > next_element_level or l not in dir_levels[n + 1:]:
                            
                            new_element = '\u2514' + '\u2500' * (n_spaces - 2) + ' ' + basename
                        
                        else:
                            
                            new_element = '\u251c' + '\u2500' * (n_spaces - 2) + ' ' + basename
                    
                    line_elements.append(new_element)
                
                structure.append(''.join(line_elements))
            
        structure = '\n'.join(structure)
        
        return structure
    
    ### ------------------------------------ ###
    
    def makeBody(self, lines):
        
        body = ['<body>']
        
        for line in lines:
            
            if line[:7] == 'HEADER:' or line[:10] == 'SUBHEADER:': # Adding (sub)header
                
                new_body_element = self.addHeader(line, self.header_counts)
                body.append(new_body_element)
                self.header_counts += 1
            
            elif line[:5] != 'PLOT:' and line[:6] != 'TABLE:': # Paragraph
                
                new_body_element = self.addParagraph(line)
                body.append(new_body_element)
            
            elif line[:5] == 'PLOT:': # Plot
                
                new_body_element = self.addPlot(line)
                body.append(new_body_element)
                
            elif line[:6] == 'TABLE:': # Table
                
                new_body_element = self.addTable(line)
                body.append(new_body_element)
            
            else:
                
                next
        
        body.append('</body>')
        body = '\n'.join(body)
        
        return body
    
    ### ------------------------------------ ###
    
    @staticmethod
    def addHeader(title, index):
        
        if title[:7] == 'HEADER:':
            
            header = f"""<h1 id="ch{index}" style="color: white; background-color: #4c4cff; border: 2px solid black; border-radius: 5px 5px 5px 5px;">{title.replace('HEADER:', '')}</h1></head>"""
        
        else:
            
            header = f"""<h2 id="ch{index}" style="color: white; background-color: #4c97ff; border: 2px solid black; border-radius: 5px 5px 5px 5px;">{title.replace('SUBHEADER:', '')}</h2></head>"""
        
        return header
    
    ### ------------------------------------ ###
    
    @staticmethod
    def addParagraph(text):
        
        paragraph = f"""
        <p style="text-align:justify; font-size: 20px; background-color: white;">
        {text}
        </p>
        """
        
        return paragraph
    
    ### ------------------------------------ ###
    
    @staticmethod
    def addPlot(info):
        
        # Plot title and path
        title, path = info.split(':')[1:]
        
        # Adding plots (1 per line)
        plot_field = f"""
        <div style="text-align:center; font-size: 25px;">
        <figcaption><br><b>{title}</b></figcaption>
        <img src="{path}" width="1024"/>
        </div>
        """
        
        return plot_field
    
    ### ------------------------------------ ###
    
    @staticmethod
    def addTable(info):
        
        # Plot title and path
        title, path = info.split(':')[1:]
        
        # Import table, then save it to html-formatted string
        table = pd.read_csv(path, sep='\t')
        table = f"""
        <p style="text-align:center; font-size: 25px;"><br><b>{title}</b></p>
        <div style='border:2px solid black; height: fit-content; min-height: 250px; max-height: 500px; overflow: auto; width: fit-content; min-width: 512; max-width: 1024px; margin-left: auto; margin-right: auto;'>
        {table.to_html(justify='center')}
        </div>
        <p><br></p>
        """
        
        return table

### ------------------MAIN------------------ ###

import pandas as pd

from os import chdir, walk
from sys import argv

try:

    basedir = argv[argv.index('--basedir') + 1]

except:
    
    basedir = './'

chdir(basedir)

try:

    readme = argv[argv.index('--readme_file') + 1]

except:
    
    readme = 'ReadME.txt'

_ = produceReport(readme)
