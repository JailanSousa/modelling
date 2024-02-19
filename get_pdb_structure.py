import requests

def pdb_template():
    
    url = 'https://files.rcsb.org/view/1CN2.pdb'
    headers = {
        'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
        'Accept-Encoding':'gzip, deflate, br',
        'Accept-Language':'pt-BR,pt;q=0.9,en-US;q=0.8,en;q=0.7',
        'User-Agent':"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
        }

    request = requests.get(url=url, headers=headers)
    return request.text
