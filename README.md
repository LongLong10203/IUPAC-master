# IUPAC Master

A website that helps you learn [IUPAC naming](https://en.wikipedia.org/wiki/IUPAC_nomenclature_of_organic_chemistry), designed for everyone.

## How to Play
1. Click on the [link](https://iupac-master.xagentx.link).
2. Enter your name and password, the website will automatically log you in.
3. Click the `NEW GAME` button.
4. A random organic compound will be generated. Enter the IUPAC name of the compound. Submit by pressing `Enter`.
5. Click `Next` to continue. If you get it right, your score will be increased. If you get it wrong, you will be redirected to the result page.
6. Compete with your friends and see who can get the highest score!

## Important Note
Please do NOT attempt to cheat. Cheating deprives you of the opportunity to learn and understand IUPAC naming, which is the goal of this game. Let's play fair and have fun!

## Libraries and Frameworks Used
- **[Flask](https://flask.palletsprojects.com/):** A lightweight WSGI web application framework in Python, used for developing web applications.
- **[Jinja](https://jinja.palletsprojects.com/):** A template engine for Python, used in Flask for rendering dynamic HTML templates.
- **[RDKit](https://www.rdkit.org/):** A collection of cheminformatics tools for molecular modeling, representation, and computation.
- **[PubChemPy](https://pubchempy.readthedocs.io/en/latest/):** A Python library for interacting with the PubChem chemical database
- **[Prisma](https://www.prisma.io/):** An experimental Object Relational Mapping (ORM) system with Python support, used for managing database interactions.
- **[PostgreSQL](https://www.postgresql.org/):** A robust, open-source object-relational database system, compatible with Prisma for database management.

## Special Thanks to
- Ms. Cheung for teaching me chemistry (although I'm failing).
- My friend [Jay](https://github.com/Agent-01) for hosting the website with his server.
- [Coolify](https://coolify.io/) for the cloud service.
- [UIverse](https://uiverse.io/) for the wonderful UI designs.

## License
[MIT License](https://en.wikipedia.org/wiki/MIT_License)

## Update Patch
- 24/11/2024: Demo release
- 25/11/2024: Alpha v1 version release
    - Secure data transmission
    - Add scoreboard page
    - Set up database to store users' password and maximum scores
    - Improve UI designs
- 26/11/2024: Alpha v2 version release
    - Bug fix
- 27/11/2024: Beta v1 release
    - Bug fix
    - Improve user experience
    - **Redesign the algorithm for generating organic compounds**
- 29/11/2024: Beta v2 release
    - Vulnerability fix to prevent cheating
    - Update database schema to include created_at and updated_at column
    - Add "Save & Exit" button to the game page

### The website has not been officially released yet