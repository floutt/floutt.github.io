AUTHOR = 'Oluwatosin Olayinka'
HANDLE = 'floutt'
SITENAME = 'floutt'
SITEURL = ''                      # empty for local dev; publishconf.py sets the real URL
GITHUB_URL = 'https://github.com/floutt'   # <- change if your username differs

PATH = 'content'
TIMEZONE = 'America/New_York'
DEFAULT_LANG = 'en'

THEME = 'theme'
THEME_STATIC_DIR = 'theme'

# Home = /, blog list = /blog/, about = /about/
DIRECT_TEMPLATES = ['index', 'blog']
BLOG_SAVE_AS = 'blog/index.html'
BLOG_URL = 'blog/'

ARTICLE_URL = 'blog/{slug}/'
ARTICLE_SAVE_AS = 'blog/{slug}/index.html'
PAGE_URL = '{slug}/'
PAGE_SAVE_AS = '{slug}/index.html'

# Only generate what the site uses
ARCHIVES_SAVE_AS = ''
AUTHORS_SAVE_AS = ''
AUTHOR_SAVE_AS = ''
CATEGORIES_SAVE_AS = ''
CATEGORY_SAVE_AS = ''
TAGS_SAVE_AS = ''
TAG_SAVE_AS = ''
PERIOD_ARCHIVE_SAVE_AS = ''

# Feeds (Atom)
FEED_DOMAIN = SITEURL
FEED_ALL_ATOM = 'feeds/atom.xml'
FEED_ALL_RSS = None
CATEGORY_FEED_ATOM = None
TRANSLATION_FEED_ATOM = None
AUTHOR_FEED_ATOM = None
AUTHOR_FEED_RSS = None

DEFAULT_PAGINATION = False
STATIC_PATHS = ['images', 'extra']
EXTRA_PATH_METADATA = {
    'extra/CNAME': {'path': 'CNAME'},
    'extra/.nojekyll': {'path': '.nojekyll'},
}
RELATIVE_URLS = False
