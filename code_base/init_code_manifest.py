#!/usr/bin/env python3

### ---------------------------------------- ###

### ------------------MAIN------------------ ###

import json

from os import listdir
from os.path import isdir

### Init manifest

manifest = {'scope_of_work' : [],
            'code' : {}}

### Add scripts to manifest

accepted_script_types = ['.sh', '.R', '.py']

for main_tag in listdir():
    
    if not isdir(main_tag):
        
        continue
    
    if main_tag not in manifest['scope_of_work']:
        
        manifest['scope_of_work'].append(main_tag)
    
    for secondary_tag in listdir(main_tag):
        
        if not isdir(f'{main_tag}/{secondary_tag}'):
            
            continue
        
        for code_name in listdir(f'{main_tag}/{secondary_tag}'):
            
            if sum([1 for ast in accepted_script_types if code_name.endswith(ast)]) > 0:
                
                path = f'{main_tag}/{secondary_tag}/{code_name}'
                
                structured_code = {'path' : path,
                                   'main_tag' : main_tag,
                                   'secondary_tag' : secondary_tag}
                
                manifest['code'][code_name] = structured_code

### Add github code

if 'github_code.json' in listdir():

    with open('github_code.json', 'r') as github_code_open:
    
        github_code = json.load(github_code_open)
    
    for code_name,code_info in github_code.items():
        
        if code_info['main_tag'] not in manifest['scope_of_work']:
            
            manifest['scope_of_work'].append(code_info['main_tag'])
        
        manifest['code'][code_name] = code_info

### Sort based on scope of work

manifest['scope_of_work'].sort()

sorting_parameters = [[i['main_tag'], i['secondary_tag'], n] for n,i in manifest['code'].items()]
sorting_parameters.sort()
manifest['code'] = {code_name : manifest['code'][code_name] for (_,_,code_name) in sorting_parameters}

### Export manifest

with open('code_manifest.json', 'w') as manifest_out:
        
    json.dump(manifest, manifest_out, indent=4)
