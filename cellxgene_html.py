import re
import scanpy as sc
import json

        
def location_block(key, port):
    return """
        location /cellxgene/%s/ {
            proxy_set_header X-Forwarded-Host $host:$server_port;
            proxy_set_header X-Forwarded-Proto $scheme;
            proxy_set_header X-Forwarded-Server $host;
            proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
            proxy_pass http://localhost:%d/;
        }
    """% (key, port)

# print(location_block('macro', 5006))

config_file = '/usr/local/etc/nginx/nginx.conf'
insert_index = 0
keys_existed = []
with open(config_file, 'r') as inh:
    lines = [line.rstrip() for line in inh]
    for i, line in enumerate(lines):
        if re.search('/Users/feifei/Projects/cellxgene', line):
            insert_index = i + 3
        
        if re.search('location /cellxgene/.*?/ {', line):
            key = re.search('location /cellxgene/(.*?)/ {', line).groups()[0]
            keys_existed.append(key)

print(keys_existed)

# read info file from json to dict
cellxgene_data = None
with open('cellxgene.json', 'r') as inh:
    cellxgene_data=json.load(inh)

    # for key, v in data.items():
    #     title = v['title']
    #     tissue = v['tissue']
    #     assay = v['assay']
    #     disease = v['disease']
    #     organism = v['organism']

            
items = []
location_blocks = ''    
with open('cellxgene.sh') as inh:
    for line in inh:
        if line.startswith('cellxgene launch'):
            arr = line.split()
            port = int(arr[5])
            file = arr[8]
            adata = sc.read_h5ad(file)
            n_cells = adata.n_obs
            key = file.split('.')[0].split('/')[1]
            title = cellxgene_data[key]['title']
            tissue = cellxgene_data[key]['tissue']
            assay = cellxgene_data[key]['assay']
            disease = cellxgene_data[key]['disease']
            organism = cellxgene_data[key]['organism']
            item = dict(key=key, title=title, tissue=tissue, assay=assay, disease=disease, organism=organism, n_cells=n_cells)
            items.append(item)
            
            # location block, pass if keys existed already
            if key not in keys_existed:
                location_blocks += location_block(key, port)

# HTML from template

from jinja2 import Environment, FileSystemLoader

file_loader = FileSystemLoader('templates')
env = Environment(loader=file_loader)

template = env.get_template('cellxgene.html')

html_output = template.render(items=items)
with open('index.html', 'w') as outh:
    outh.write(html_output)


# Update config file
new_lines = []
any_update = False
with open(config_file, 'r') as inh:
    lines = [line.rstrip() for line in inh]
    new_lines = lines
    
    # If the block key exists, then it should pass ...
    if location_blocks:
        any_update = True
        new_lines.insert(insert_index, location_blocks)
    

if any_update:
    config_file_new = 'nginx.conf'
    with open(config_file_new, 'w') as outh:
        outh.write('\n'.join(new_lines))
