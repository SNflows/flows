from tendrils import utils
import os


def set_up_tendrils_in_ci():
    # Set UP API key
    os.environ.get('FLOWS_API_TOKEN')
    utils.set_api_token(token=os.environ.get('FLOWS_API_TOKEN'), overwrite=True)


if __name__ == '__main__':
    set_up_tendrils_in_ci()
