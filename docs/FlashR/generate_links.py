import glob
import ntpath

DIRECTORY = "."
ROOT_URL = "http://flashx.io/docs/FlashR/"

links = []
for f in glob.glob('{}/*.html'.format(DIRECTORY)):
    f = ntpath.basename(f)
    links.append("<li><a href='{}{}'>{}</a></li>".format(ROOT_URL, f, f))

print("<ul>" + "\n".join(links) + "</ul>")
