import logging

__version__ = "1.4.2"

# Library convention: attach a NullHandler so trviz never emits log records
# unless the user explicitly configures logging in their application. See
# https://docs.python.org/3/howto/logging.html#configuring-logging-for-a-library
logging.getLogger(__name__).addHandler(logging.NullHandler())
