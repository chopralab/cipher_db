class UserNotFoundError(Exception):
    pass

class LoginFailedError(Exception):
    pass

class InvalidUsernameError(Exception):
    pass

class InvalidPasswordError(Exception):
    pass

def login(username, password):
    '''
    Attemps to log a user into the website

    Parameters
    ----------
    username: string, required
        The username of the specified user
    password: string, required
        The password of the specified user

    Raises
    ------
    LoginFailedError
        If the login credentials do not match
    '''
    pass

def add_user(username, password):
    '''
    Adds a new user to the database with the specified username and password

    Parameters
    ----------
    username: string, required
        The username of the user to be added
    password: string, required
        The password of the user to be added

    Raises
    ------
    InvalidUsernameError
        If the username does not fit the creation criteria
    InvalidPasswordError
        If the password does not fit the creation criteria
    '''
    pass

def remove_user(username, password):
    '''
    Removes the specified user from the database

    Parameters
    ----------
    username: string, required
        The username of the specified user
    password: string, required
        The password of the specified user
    
    Raises
    ------
    UserNotFoundError
        If the user is not found in the database
    '''
    pass

def change_username(username, password, new_username):
    '''
    Changes the username of the specified user

    Parameters
    ----------
    username: string, required
        The current username of the specified user
    password: string, required
        The password of the specified user
    new_username: string, required
        The new username of the specified user

    Raises
    ------
    UserNotFoundError
        If the user is not found in the database
    LoginFailedError
        If the login credentials do not match
    InvalidUsernameError
        If the username does not fit the creation criteria
    '''
    pass

def change_password(username, password, new_password):
    '''
    Changes the password of the specified user

    Parameters
    ----------
    username: string, required
        The username of the specified user
    password: string, required
        The current password of the specified user
    new_password: string, required
        The new password of the specified user

    Raises
    ------
    UserNotFoundError
        If the user is not found in the database
    LoginFailedError
        If the login credentials do not match
    InvalidPasswordError
        If the password does not fit the creation criteria
    '''
    pass