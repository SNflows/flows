from tendrils import utils
import os


def set_up_tendrils_in_ci():
    # Set UP API key
    config = utils.load_config()
    if config.has_option('api', 'token'):
        token = config.get('api', 'token')
        if token not in [None, '', 'None', 'test']:  # if not, we probably have a token.
            return

    os.environ.get('FLOWS_API_TOKEN')
    utils.set_api_token(token=os.environ.get('FLOWS_API_TOKEN'), overwrite=True)


if __name__ == '__main__':
    set_up_tendrils_in_ci()
