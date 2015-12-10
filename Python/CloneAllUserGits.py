# coding: utf-8

"""
modified from https://filippo.io/archive-your-github-repo-and-data/
to allow up to 'n' repos ('n_repos'). Git API defaults out at 30 repos.
"""

import argparse
from urllib import urlopen
from subprocess import call
import json
import re
import os.path

parser = argparse.ArgumentParser(description='Dump an user\'s public GitHub data into current directory.')
parser.add_argument('-u', '--user', required=True, help='username')
parser.add_argument('-n', '--n_repos', default = 30, help = 'Number of repos to stop downloading at (Git API defaults out at 30 – if you have a large number make this something like 100 or 1000 or even 5000 to be sure).', required=False)

args = parser.parse_args()

def clear_url(url):
    return re.sub(r'\{[^\}]*\}', '', url)

data = urlopen('https://api.github.com/users/' + args.user).read()
user = json.loads(data.decode('utf-8'))
data = urlopen(clear_url(user['repos_url']+'?per_page='+str(args.n_repos))).read()
repos = json.loads(data.decode('utf-8'))
for repo in repos:
    if not repo['fork']:
        call(['git', 'clone', repo['clone_url']])
    elif args.forks:
        if not os.path.exists('forks'):
            os.makedirs('forks')
        call(['git', 'clone', repo['clone_url'], os.path.join('forks', repo['name'])])
