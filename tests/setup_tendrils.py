from tendrils import utils
import os


def set_up_tendrils_in_ci():
    # Set UP API key
    config = utils.load_config()
    if config.has_option('api', 'token'):
        token = config.get('api', 'token')
        if token not in [None, '', 'None', 'test']:  # if not, we probably have a token.
            return

    token = os.environ.get('FLOWS_API_KEY')
    if token is None:
        raise RuntimeError("API token can not be set up.")
    utils.set_api_token(token=token, overwrite=True)


if __name__ == '__main__':
    set_up_tendrils_in_ci()
